#include "gpuRetroGatedSensePrepGadget.h"
#include "cuNonCartesianSenseOperator.h"
#include "GenericReconJob.h"
#include "cuNDArray_elemwise.h"
#include "hoNDArray_elemwise.h"
#include "cuNDArray_utils.h"
#include "vector_td_operators.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "check_CUDA.h"
#include "radial_utilities.h"
#include "hoNDArray_fileio.h"
#include "ismrmrd/xml.h"

#include <algorithm>
#include <vector>
#include <cmath>

namespace Gadgetron{

  gpuRetroGatedSensePrepGadget::gpuRetroGatedSensePrepGadget()
    : slices_(-1)
    , sets_(-1)
    , samples_per_profile_(-1)
    , phys_time_index_(0)
  {

  }
  
  gpuRetroGatedSensePrepGadget::~gpuRetroGatedSensePrepGadget() {}
  
  int gpuRetroGatedSensePrepGadget::process_config(ACE_Message_Block* mb)
  {
    // Get configuration values from config file
    //

    mode_ = mode.value();
    device_number_ = deviceno.value();
    profiles_per_frame_ = profiles_per_frame.value();
    frames_per_cardiac_cycle_ = frames_per_cardiac_cycle.value();
    profiles_per_buffer_frame_ = profiles_per_buffer_frame.value();
    num_buffer_frames_inner_ = number_of_buffer_frames_inner.value();
    num_buffer_frames_outer_ = number_of_buffer_frames_outer.value();
    buffer_using_solver_ = buffer_using_solver.value();
    output_timing_ = output_timing.value();
    phys_time_index_ = physiology_time_index.value();

    // Check that a golden ratio based reconstruction mode was specified
    //

    if( !(mode_ == 2 || mode_ == 3) ){
      GDEBUG( "Only radial reconstruction modes {2,3} (golden ratio based) are supported.\n" );
      return GADGET_FAIL;
    }
    
    // Setup and validate device configuration
    //

    int number_of_devices;
    if (cudaGetDeviceCount(&number_of_devices)!= cudaSuccess) {
      GDEBUG( "Error: unable to query number of CUDA devices.\n" );
      return GADGET_FAIL;
    }

    if (number_of_devices == 0) {
      GDEBUG( "Error: No available CUDA devices.\n" );
      return GADGET_FAIL;
    }

    if (device_number_ >= number_of_devices) {
      GDEBUG("Adjusting device number from %d to %d\n", device_number_,  (device_number_%number_of_devices));
      device_number_ = (device_number_%number_of_devices);
    }

    if (cudaSetDevice(device_number_)!= cudaSuccess) {
      GDEBUG( "Error: unable to set CUDA device.\n" );
      return GADGET_FAIL;
    }

    cudaDeviceProp deviceProp;
    if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
      GDEBUG( "Error: unable to query device properties.\n" );
      return GADGET_FAIL;
    }
    
    // Convolution kernel width and oversampling ratio (for the buffer)
    //

    kernel_width_ = buffer_convolution_kernel_width.value();
    oversampling_factor_ = buffer_convolution_oversampling_factor.value();

    // Get the Ismrmrd header
    //

    ISMRMRD::IsmrmrdHeader h;
    ISMRMRD::deserialize(mb->rd_ptr(),h);
    
    
    if (h.encoding.size() != 1) {
      GDEBUG("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }
    
    // Get the encoding space and trajectory description
    ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
    ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
    ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;


    // Matrix sizes (as a multiple of the GPU's warp size)
    //
    
    unsigned int warp_size = deviceProp.warpSize;

    image_dimensions_.push_back(((e_space.matrixSize.x+warp_size-1)/warp_size)*warp_size);
    image_dimensions_.push_back(((e_space.matrixSize.y+warp_size-1)/warp_size)*warp_size);

    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(e_space.matrixSize.x*reconstruction_os_factor_x.value()))+warp_size-1)/warp_size)*warp_size);  
    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(e_space.matrixSize.y*reconstruction_os_factor_y.value()))+warp_size-1)/warp_size)*warp_size);
    
    image_dimensions_recon_os_ = uint64d2
      (((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
       ((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);
    
    // In case the warp_size constraint kicked in
    oversampling_factor_ = float(image_dimensions_recon_os_[0])/float(image_dimensions_recon_[0]); 
    
    GDEBUG("matrix_size_x : %d, recon: %d, recon_os: %d\n", 
                  image_dimensions_[0], image_dimensions_recon_[0], image_dimensions_recon_os_[0]);

    GDEBUG("matrix_size_y : %d, recon: %d, recon_os: %d\n", 
                  image_dimensions_[1], image_dimensions_recon_[1], image_dimensions_recon_os_[1]);
    
    fov_.push_back(r_space.fieldOfView_mm.x);
    fov_.push_back(r_space.fieldOfView_mm.y);
    fov_.push_back(r_space.fieldOfView_mm.z);

    slices_ = e_limits.slice ? e_limits.slice->maximum + 1 : 1;
    sets_ = e_limits.set ? e_limits.set->maximum + 1 : 1;
    
    // Allocate profile queues
    // - one queue for the currently incoming frame (for the accumulation buffer)
    // - one queue for the next reconstruction
    
    buffer_profiles_queue_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);
    recon_profiles_queue_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);

    size_t bsize = sizeof(GadgetContainerMessage< hoNDArray< std::complex<float> > >)*profiles_per_buffer_frame_*10;

    for( unsigned int i=0; i<slices_*sets_; i++ ){
      buffer_profiles_queue_[i].high_water_mark(bsize);
      buffer_profiles_queue_[i].low_water_mark(bsize);
    }

    bsize = sizeof(GadgetContainerMessage< hoNDArray< std::complex<float> > >)*profiles_per_frame_*frames_per_cardiac_cycle_*10;
    
    for( unsigned int i=0; i<slices_*sets_; i++ ){
      recon_profiles_queue_[i].high_water_mark(bsize);
      recon_profiles_queue_[i].low_water_mark(bsize);
    }
    
    // Define some profile counters for book-keeping
    //

    image_counter_ = boost::shared_array<long>(new long[slices_*sets_]);
    num_coils_ = boost::shared_array<unsigned int>(new unsigned int[slices_*sets_]);
    first_profile_acq_time_ = boost::shared_array<unsigned int>(new unsigned int[slices_*sets_]);
    first_profile_phys_time_ = boost::shared_array<unsigned int>(new unsigned int[slices_*sets_]);
    previous_timestamp_ = boost::shared_array<unsigned int>(new unsigned int[slices_*sets_]);
    profiles_counter_global_ = boost::shared_array<long>(new long[slices_*sets_]);
    Rw_reached_ = boost::shared_array<bool>(new bool[slices_*sets_]);
    Rw_offset_ = boost::shared_array<unsigned int>(new unsigned int[slices_*sets_]);
    buffer_update_needed_ = boost::shared_array<bool>(new bool[slices_*sets_]);
    reconfigure_ = boost::shared_array<bool>(new bool[slices_*sets_]);
    
    if( !image_counter_.get() || 
        !num_coils_.get() || 
        !first_profile_acq_time_.get() ||
        !first_profile_phys_time_.get() ||
        !previous_timestamp_.get() ||
        !profiles_counter_global_.get() ||
        !Rw_reached_.get() ||
        !Rw_offset_.get() ||
        !buffer_update_needed_.get() ||
        !reconfigure_ ){
      GDEBUG("Failed to allocate host memory (1)\n");
      return GADGET_FAIL;
    }

    for( unsigned int i=0; i<slices_*sets_; i++ ){
      image_counter_[i] = 0;
      num_coils_[i] = 0;
      previous_timestamp_[i] = 0;
      profiles_counter_global_[i] = 0;
      Rw_reached_[i] = false;
      Rw_offset_[i] = 0;
      buffer_update_needed_[i] = true;
      reconfigure_[i] = true;
    }
        
    position_ = boost::shared_array<float[3]>(new float[slices_*sets_][3]);
    read_dir_ = boost::shared_array<float[3]>(new float[slices_*sets_][3]);
    phase_dir_ = boost::shared_array<float[3]>(new float[slices_*sets_][3]);
    slice_dir_ = boost::shared_array<float[3]>(new float[slices_*sets_][3]);

    if( !position_.get() || !read_dir_.get() || !phase_dir_.get() || !slice_dir_.get() ){
      GDEBUG("Failed to allocate host memory (2)\n");
      return GADGET_FAIL;
    }

    for( unsigned int i=0; i<slices_*sets_; i++ ){
      (position_[i])[0] = (position_[i])[1] = (position_[i])[2] = 0.0f;
      (read_dir_[i])[0] = (read_dir_[i])[1] = (read_dir_[i])[2] = 0.0f;
      (phase_dir_[i])[0] = (phase_dir_[i])[1] = (phase_dir_[i])[2] = 0.0f;
      (slice_dir_[i])[0] = (slice_dir_[i])[1] = (slice_dir_[i])[2] = 0.0f;
    }

    // Allocate accumulation buffer
    //

    if( buffer_using_solver_ )
      acc_buffer_cg_ = boost::shared_array< cuSenseBufferCg<float,2> >(new cuSenseBufferCg<float,2>[slices_*sets_]);
    else
      acc_buffer_ = boost::shared_array< cuSenseBuffer<float,2> >(new cuSenseBuffer<float,2>[slices_*sets_]);
    
    // Allocate remaining shared_arrays
    //
    
    csm_host_ = boost::shared_array< hoNDArray<float_complext> >(new hoNDArray<float_complext>[slices_*sets_]);
    reg_host_ = boost::shared_array< hoNDArray<float_complext> >(new hoNDArray<float_complext>[slices_*sets_]);

    host_weights_recon_ = boost::shared_array< hoNDArray<float> >(new hoNDArray<float>[slices_*sets_]);

    if( !csm_host_.get() || !reg_host_.get() || !host_weights_recon_ ){
      GDEBUG("Failed to allocate host memory (3)\n");
      return GADGET_FAIL;
    }

    return GADGET_OK;
  }

  int gpuRetroGatedSensePrepGadget::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
          GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2)
  {
    // Noise should have been consumed by the noise adjust (if in the gadget chain)
    //
    
    bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
    if (is_noise) { 
      m1->release();
      return GADGET_OK;
    }

    unsigned int slice = m1->getObjectPtr()->idx.slice;
    unsigned int set = m1->getObjectPtr()->idx.set;

    unsigned int profile = m1->getObjectPtr()->idx.kspace_encode_step_1;

    unsigned int current_timestamp = m1->getObjectPtr()->physiology_time_stamp[phys_time_index_];
    unsigned int previous_timestamp = previous_timestamp_[set*slices_+slice];
    
    bool new_cardiac_cycle_detected = (current_timestamp < previous_timestamp);

    previous_timestamp_[set*slices_+slice] = current_timestamp;

    if( !Rw_reached_[set*slices_+slice] && !new_cardiac_cycle_detected ){ 
      Rw_offset_[set*slices_+slice]++;
      m1->release();
      return GADGET_OK;
    }

    if( !Rw_reached_[set*slices_+slice] && new_cardiac_cycle_detected ){ 
      Rw_reached_[set*slices_+slice] = true;
      profiles_counter_global_[set*slices_+slice] = Rw_offset_[set*slices_+slice];
      new_cardiac_cycle_detected = false;
    }

    boost::shared_ptr<GPUTimer> process_timer;
    if( output_timing_ )
      process_timer = boost::shared_ptr<GPUTimer>( new GPUTimer("gpuRetroGatedSensePrepGadget::process()") );

    // Get a pointer to the accumulation buffer. 
    //

    cuSenseBuffer<float,2> *acc_buffer = (buffer_using_solver_) ? &acc_buffer_cg_[set*slices_+slice] : &acc_buffer_[set*slices_+slice];

    // Have the imaging plane changed?
    //

    if( !vec_equal(position_[set*slices_+slice], m1->getObjectPtr()->position) ||
        !vec_equal(read_dir_[set*slices_+slice], m1->getObjectPtr()->read_dir) || 
        !vec_equal(phase_dir_[set*slices_+slice], m1->getObjectPtr()->phase_dir) ||
        !vec_equal(slice_dir_[set*slices_+slice], m1->getObjectPtr()->slice_dir) ){
      
      // Yes indeed, clear the accumulation buffer
      acc_buffer->clear();
      buffer_update_needed_[set*slices_+slice] = true;
      
      memcpy(position_[set*slices_+slice],m1->getObjectPtr()->position,3*sizeof(float));
      memcpy(read_dir_[set*slices_+slice],m1->getObjectPtr()->read_dir,3*sizeof(float));
      memcpy(phase_dir_[set*slices_+slice],m1->getObjectPtr()->phase_dir,3*sizeof(float));
      memcpy(slice_dir_[set*slices_+slice],m1->getObjectPtr()->slice_dir,3*sizeof(float));
    }
    
    // Only when the first profile arrives, do we know the #samples/profile
    //

    if( samples_per_profile_ == -1 )      
      samples_per_profile_ = m1->getObjectPtr()->number_of_samples;
    
    if( samples_per_profile_ != m1->getObjectPtr()->number_of_samples ){
      GDEBUG("Unexpected change in the incoming profiles' lengths\n");
      return GADGET_FAIL;
    }
    
    // Reconfigure at first pass
    // - or if the number of coil changes
    // - or if the reconfigure_ flag is set

    if( num_coils_[set*slices_+slice] != m1->getObjectPtr()->active_channels ){
      GDEBUG("Reconfiguring due to change in the number of coils\n");
      num_coils_[set*slices_+slice] = m1->getObjectPtr()->active_channels;
      reconfigure(set, slice);
    }

    if( reconfigure_[set*slices_+slice] ){
      GDEBUG("Reconfiguring due to boolean indicator\n");
      reconfigure(set, slice);
    }

    // Enqueue profile
    // - if 'new_cardiac_cycle_detected' the current profile does not
    //   belong to the current cardiac cycle and we delay enqueing
    //

    buffer_profiles_queue_[set*slices_+slice].enqueue_tail(duplicate_profile(m2));
    
    if( !new_cardiac_cycle_detected ) {
      if( recon_profiles_queue_[set*slices_+slice].message_count() == 0 ){
        first_profile_acq_time_[set*slices_+slice] = m1->getObjectPtr()->acquisition_time_stamp;
        first_profile_phys_time_[set*slices_+slice] = m1->getObjectPtr()->physiology_time_stamp[phys_time_index_];
      }
      recon_profiles_queue_[set*slices_+slice].enqueue_tail(duplicate_profile(m2));
    }
    
    // If the profile is the last of a "buffer frame" 
    // - then update the accumulation buffer
    //
    
    bool is_last_profile_in_buffer_frame = 
      ( buffer_profiles_queue_[set*slices_+slice].message_count() == profiles_per_buffer_frame_ );
    
    if( is_last_profile_in_buffer_frame ){
      
      // Extract this frame's samples to update the csm/regularization buffer
      //
      
      boost::shared_ptr< hoNDArray<float_complext> > host_samples = 
        extract_samples_from_buffer_queue( set, slice );
      
      if( host_samples.get() == 0x0 ){
        GDEBUG("Failed to extract buffer samples from queue\n");
        return GADGET_FAIL;
      }
      
      cuNDArray<float_complext> samples( host_samples.get() );
      
      long profile_offset = profiles_counter_global_[set*slices_+slice];
      boost::shared_ptr< cuNDArray<floatd2> > traj = calculate_trajectory_for_buffer(profile_offset, set, slice);
      
      buffer_update_needed_[set*slices_+slice] |= acc_buffer->add_frame_data( &samples, traj.get() );
    }
    
    // Perform reconstruction if it is time...
    //
      
    if( new_cardiac_cycle_detected ){
      
      // Prepare the image headers for the reconstruction
      //
      
      boost::shared_array<ISMRMRD::ImageHeader> headers( new ISMRMRD::ImageHeader[frames_per_cardiac_cycle_] );
      
      for( unsigned int i=0; i<frames_per_cardiac_cycle_; i++ ){
        
        ISMRMRD::AcquisitionHeader *base_head = m1->getObjectPtr();
        ISMRMRD::ImageHeader *header = &headers[i];
        
        {
          // Initialize header to all zeroes (there is a few fields we do not set yet)
          ISMRMRD::ImageHeader tmp;
          *header = tmp;
        }
        
        header->version = base_head->version;
        
        header->matrix_size[0] = image_dimensions_recon_[0];
        header->matrix_size[1] = image_dimensions_recon_[1];
        header->matrix_size[2] = 1;
        
        header->field_of_view[0] = fov_[0];
        header->field_of_view[1] = fov_[1];
        header->field_of_view[2] = fov_[2];
        
        header->channels = num_coils_[set*slices_+slice];
        header->slice = base_head->idx.slice;
        header->set = base_head->idx.set;
        
        header->acquisition_time_stamp = 
          first_profile_acq_time_[set*slices_+slice] + 
          i*(base_head->acquisition_time_stamp-first_profile_acq_time_[set*slices_+slice])/frames_per_cardiac_cycle_;

        header->physiology_time_stamp[phys_time_index_] = 
          first_profile_phys_time_[set*slices_+slice] + 
          i*(base_head->physiology_time_stamp[phys_time_index_]-first_profile_phys_time_[set*slices_+slice])/frames_per_cardiac_cycle_;

        memcpy(header->position, base_head->position, sizeof(float)*3);
        memcpy(header->read_dir, base_head->read_dir, sizeof(float)*3);
        memcpy(header->phase_dir, base_head->phase_dir, sizeof(float)*3);
        memcpy(header->slice_dir, base_head->slice_dir, sizeof(float)*3);
        memcpy(header->patient_table_position, base_head->patient_table_position, sizeof(float)*3);
        
        header->data_type = ISMRMRD::ISMRMRD_CXFLOAT;
        header->image_index = image_counter_[set*slices_+slice]++; 
        header->image_series_index = set*slices_+slice;        
      }
      
      // Update csm and regularization images
      //

      if( buffer_update_needed_[set*slices_+slice] || 
          csm_host_[set*slices_+slice].get_number_of_elements() == 0 || 
          reg_host_[set*slices_+slice].get_number_of_elements() == 0 ) {

        // Get the accumulated coil images
        //
        
        boost::shared_ptr< cuNDArray<float_complext> > csm_data = acc_buffer->get_accumulated_coil_images();
        
        if( !csm_data.get() ){
          GDEBUG("Error during accumulation buffer computation\n");
          return GADGET_FAIL;
        }            
	
        // Estimate CSM
        //

        auto csm = boost::make_shared<cuNDArray<float_complext>>(estimate_b1_map<float,2>( csm_data.get() ));


        acc_buffer->set_csm(csm);
        csm_host_[set*slices_+slice] = *(csm->to_host());
      
        // Compute regularization image
        //

        boost::shared_ptr< cuNDArray<float_complext> > reg_image;
	
        if( buffer_using_solver_ ){
          ((cuSenseBufferCg<float,2>*)acc_buffer)->preprocess( calculate_trajectory_for_rhs( profiles_counter_global_[set*slices_+slice], set, slice).get() );
        }
      
        reg_image = acc_buffer->get_combined_coil_image();
        
        if( !reg_image.get() ){
          GDEBUG("Error computing regularization image\n");
          return GADGET_FAIL;
        }            
	
        reg_host_[set*slices_+slice] = *(reg_image->to_host());
        
        /*
          static int counter = 0;
          char filename[256];
          sprintf((char*)filename, "reg_%d.real", counter);
          write_nd_array<float>( abs(&reg_host_[set*slices_+slice]).get(), filename );
          counter++;  */

        buffer_update_needed_[set*slices_+slice] = false;        
      }

      // Prepare data array of the profiles for the downstream reconstruction
      //
      
      boost::shared_ptr< hoNDArray<float_complext> > samples_host( new hoNDArray<float_complext>() );
      boost::shared_ptr< hoNDArray<floatd2> > traj_host( new hoNDArray<floatd2> );

      if( extract_samples_and_trajectory_from_recon_queue( set, slice, samples_host, traj_host ) != GADGET_OK ){
        GDEBUG("Failed to extract samples and/or trajectories.\n");
        return GADGET_FAIL;
      }        
      
      // Set up Sense job
      //

      GadgetContainerMessage< GenericReconJob >* m4 = new GadgetContainerMessage< GenericReconJob >();
	
      m4->getObjectPtr()->dat_host_ = samples_host;
      m4->getObjectPtr()->tra_host_ = traj_host;
      m4->getObjectPtr()->dcw_host_ = boost::shared_ptr< hoNDArray<float> >(new hoNDArray<float>(host_weights_recon_[set*slices_+slice]));
      m4->getObjectPtr()->csm_host_ = boost::shared_ptr< hoNDArray<float_complext> >( new hoNDArray<float_complext>(csm_host_[set*slices_+slice]));
      m4->getObjectPtr()->reg_host_ = boost::shared_ptr< hoNDArray<float_complext> >( new hoNDArray<float_complext>(reg_host_[set*slices_+slice]));
      m4->getObjectPtr()->image_headers_ = headers;
      
      // The Sense Job needs an image header as well. 
      // Let us just copy the initial one...
      //

      GadgetContainerMessage<ISMRMRD::ImageHeader> *m3 = new GadgetContainerMessage<ISMRMRD::ImageHeader>;
      *m3->getObjectPtr() = m4->getObjectPtr()->image_headers_[0];
      m3->cont(m4);
      
      if (this->next()->putq(m3) < 0) {
        GDEBUG("Failed to put job on queue.\n");
        m3->release();
        return GADGET_FAIL;
      }
    }
    
    // This is was first profile of a new cardiac cycle, enqueue (since this was postponed above).
    //

    if( new_cardiac_cycle_detected ){      
      if( recon_profiles_queue_[set*slices_+slice].message_count() == 0 ){
        first_profile_acq_time_[set*slices_+slice] = m1->getObjectPtr()->acquisition_time_stamp;
        first_profile_phys_time_[set*slices_+slice] = m1->getObjectPtr()->physiology_time_stamp[phys_time_index_];
      }
      recon_profiles_queue_[set*slices_+slice].enqueue_tail(duplicate_profile(m2)); 
    }
    
    profiles_counter_global_[set*slices_+slice]++;

    if( output_timing_ )
      process_timer.reset();
    
    m1->release(); // the internal queues hold copies
    return GADGET_OK;
  }
  
  int
  gpuRetroGatedSensePrepGadget::calculate_density_compensation_for_reconstruction( unsigned int set, unsigned int slice )
  {
    switch(mode_){
      
    case 2:
    case 3:
      host_weights_recon_[set*slices_+slice] = *compute_radial_dcw_golden_ratio_2d<float>
        ( samples_per_profile_, profiles_per_frame_, oversampling_factor_, 
          1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])), 0,
          (mode_==2) ? GR_ORIGINAL : GR_SMALLEST )->to_host();
      break;
      
    default:
      GDEBUG("Illegal dcw mode\n");
      return GADGET_FAIL;
      break;
    }
    return GADGET_OK;
  }
  
  boost::shared_ptr< cuNDArray<floatd2> > 
  gpuRetroGatedSensePrepGadget::calculate_trajectory_for_buffer( long profile_offset, unsigned int set, unsigned int slice )
  {
    boost::shared_ptr< cuNDArray<floatd2> > result;

    switch(mode_){

    case 2:
    case 3:
      { 

        long first_profile_in_buffer = profile_offset + 1 - profiles_per_buffer_frame_;

        result = compute_radial_trajectory_golden_ratio_2d<float>
          ( samples_per_profile_, profiles_per_buffer_frame_, 1, first_profile_in_buffer, (mode_==2) ? GR_ORIGINAL : GR_SMALLEST );

      }
      break;	
	
    default:
      GDEBUG("Illegal trajectory mode\n");
      break;
    }
    
    return result;
  }

  boost::shared_ptr< cuNDArray<float> >
  gpuRetroGatedSensePrepGadget::calculate_density_compensation_for_buffer( unsigned int set, unsigned int slice )
  {    
    switch(mode_){
      
    case 2:
    case 3:
      return compute_radial_dcw_golden_ratio_2d<float>
        ( samples_per_profile_, profiles_per_buffer_frame_, oversampling_factor_, 
          1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])), 0,
          (mode_==2) ? GR_ORIGINAL : GR_SMALLEST );
      break;
      
    default:
      GDEBUG("Illegal dcw mode\n");
      return boost::shared_ptr< cuNDArray<float> >();
      break;
    }   
  }


  boost::shared_ptr< cuNDArray<floatd2> > 
  gpuRetroGatedSensePrepGadget::calculate_trajectory_for_rhs( long profile_offset, unsigned int set, unsigned int slice )
  {
    switch(mode_){

    case 2:
    case 3:
      { 

        long first_profile =
          std::max( 0L, profile_offset + 1 - profiles_per_buffer_frame_*num_buffer_frames_inner_ );

        return compute_radial_trajectory_golden_ratio_2d<float>
          ( samples_per_profile_, profiles_per_buffer_frame_*num_buffer_frames_inner_, 1, first_profile, (mode_==2) ? GR_ORIGINAL : GR_SMALLEST );
      }
      break;	
	
    default:
      GDEBUG("Illegal trajectory mode\n");
      return boost::shared_ptr< cuNDArray<floatd2> >();
      break;
    }
  }
  
  boost::shared_ptr< cuNDArray<float> >
  gpuRetroGatedSensePrepGadget::calculate_density_compensation_for_rhs( unsigned int set, unsigned int slice )
  {
    switch(mode_){
      
    case 2:
    case 3:
      {

        long num_profiles = profiles_per_buffer_frame_*num_buffer_frames_inner_;

        return compute_radial_dcw_golden_ratio_2d<float>
          ( samples_per_profile_, num_profiles, oversampling_factor_, 
            1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])), 0,
            (mode_==2) ? GR_ORIGINAL : GR_SMALLEST );

      }
      break;
      
    default:
      GDEBUG("Illegal dcw mode\n");
      return boost::shared_ptr< cuNDArray<float> >();
      break;
    }
  }

  boost::shared_ptr< hoNDArray<float_complext> > 
  gpuRetroGatedSensePrepGadget::extract_samples_from_buffer_queue( unsigned int set, unsigned int slice )
  {    
    ACE_Message_Queue<ACE_MT_SYNCH> *queue = &buffer_profiles_queue_[set*slices_+slice];

    unsigned int profiles_buffered = queue->message_count();
    
    std::vector<size_t> dims;
    dims.push_back(samples_per_profile_*profiles_buffered);
    dims.push_back(num_coils_[set*slices_+slice]);
    
    boost::shared_ptr< hoNDArray<float_complext> > host_samples(new hoNDArray<float_complext>(dims));
    
    for (unsigned int p=0; p<profiles_buffered; p++) {

      ACE_Message_Block* mbq;
      if (queue->dequeue_head(mbq) < 0) {
        GDEBUG("Message dequeue failed\n");
        return boost::shared_ptr< hoNDArray<float_complext> >();
      }
      
      GadgetContainerMessage< hoNDArray< std::complex<float> > > *daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq);
	
      if (!daq) {
        GDEBUG("Unable to interpret data on message queue\n");
        return boost::shared_ptr< hoNDArray<float_complext> >();
      }
	
      for (unsigned int c = 0; c < num_coils_[set*slices_+slice]; c++) {
	
        float_complext *data_ptr = host_samples->get_data_ptr();
        data_ptr += c*samples_per_profile_*profiles_buffered+p*samples_per_profile_;
	    
        std::complex<float> *r_ptr = daq->getObjectPtr()->get_data_ptr();
        r_ptr += c*daq->getObjectPtr()->get_size(0);
	  
        memcpy(data_ptr,r_ptr,samples_per_profile_*sizeof(float_complext));
      }
      
      mbq->release();
    } 
    
    return host_samples;
  }
  
  int gpuRetroGatedSensePrepGadget::extract_samples_and_trajectory_from_recon_queue
  ( unsigned int set, unsigned int slice, boost::shared_ptr< hoNDArray<float_complext> > samples, boost::shared_ptr< hoNDArray<floatd2> > trajectory )
  {    
    // Extract samples from queue and put into buffer 
    //

    ACE_Message_Queue<ACE_MT_SYNCH> *queue = &recon_profiles_queue_[set*slices_+slice];
    long profiles_buffered = queue->message_count();
    
    std::vector<size_t> dims_per_readout;
    dims_per_readout.push_back(samples_per_profile_);
    dims_per_readout.push_back(num_coils_[set*slices_+slice]);
    
    std::vector<size_t> dims_for_buffer = dims_per_readout;
    dims_for_buffer.push_back(profiles_buffered);
    
    hoNDArray< std::complex<float> > host_buffer(dims_for_buffer);

    for (long p=0; p<profiles_buffered; p++) {
      
      ACE_Message_Block* mbq;
      if (queue->dequeue_head(mbq) < 0) {
        GDEBUG("Message dequeue failed\n");
        return GADGET_FAIL;
      }
      
      GadgetContainerMessage< hoNDArray< std::complex<float> > > *daq = 
        AsContainerMessage<hoNDArray< std::complex<float> > >(mbq);
	
      if (!daq) {
        GDEBUG("Unable to interpret data on message queue\n");
        return GADGET_FAIL;
      }

      {
        // Copy daq into host_buffer array
        hoNDArray< std::complex<float> > tmp( dims_per_readout, host_buffer.get_data_ptr() + p*dims_per_readout[0]*dims_per_readout[1] );
        if( !tmp.dimensions_equal( daq->getObjectPtr()->get_dimensions().get() )){
          GDEBUG("Unexpected dimensionality of array on message queue\n");
          return GADGET_FAIL;
        }
        tmp = *daq->getObjectPtr();
      }
      mbq->release();
    } 

    // Create trajectory array according to the samples buffer
    //

    long first_profile_in_buffer = 
      profiles_counter_global_[set*slices_+slice] - profiles_buffered;
    
    boost::shared_ptr< hoNDArray<floatd2> > host_traj = compute_radial_trajectory_golden_ratio_2d<float>
      ( samples_per_profile_, profiles_buffered, 1, first_profile_in_buffer, (mode_==2) ? GR_ORIGINAL : GR_SMALLEST )->to_host();

    host_traj->squeeze();

    // Prepare samples and trajecotry arrays according to the current 
    // 'profiles_per_frame_' and 'frames_per_cardiac_cycle_' settings
    //
    
    std::vector<size_t> recon_dims;
    recon_dims.push_back(samples_per_profile_*profiles_per_frame_*frames_per_cardiac_cycle_);
    recon_dims.push_back(num_coils_[set*slices_+slice]);
    
    std::vector<size_t> traj_dims_frame;
    traj_dims_frame.push_back( samples_per_profile_*profiles_per_frame_ );

    std::vector<size_t> traj_dims = traj_dims_frame;
    traj_dims.push_back( frames_per_cardiac_cycle_ );

    samples->create( recon_dims );
    trajectory->create( traj_dims );
    
    for( long frame=0; frame<frames_per_cardiac_cycle_; frame++ ){
      
      long first_profile = 
        (long)(float(frame)*float(profiles_buffered-profiles_per_frame_)/float(frames_per_cardiac_cycle_-1));
      // Just to be sure we do run get out-of-bounds due to rounding errors in the float<->int math
      //

      if( first_profile < 0 ){
        GDEBUG("\nWARNING: first profile is negative. Corrected.");
        first_profile = 0;
      }

      if (first_profile + profiles_per_frame_ - 1  > profiles_buffered -1 ){
        GDEBUG("\nWARNING: first profile is out of bounds for the last profile. Corrected.");
        first_profile = profiles_buffered - profiles_per_frame_;
      }

      //printf( "\nFor frame %ld: The first profile has index %ld (of %ld).", frame, first_profile, profiles_buffered );
        
      for( long coil=0; coil<num_coils_[set*slices_+slice]; coil++ ){
        
        for( long profile = 0; profile<profiles_per_frame_; profile++ ){

          // Copy samples for profile
          //

          memcpy( samples->get_data_ptr() + 
                  coil*samples_per_profile_*profiles_per_frame_*frames_per_cardiac_cycle_ +
                  frame*samples_per_profile_*profiles_per_frame_ + 
                  profile*samples_per_profile_,
                  
                  host_buffer.get_data_ptr() + 
                  (first_profile + profile) * samples_per_profile_*num_coils_[set*slices_+slice]+
                  coil*samples_per_profile_,
                  
                  sizeof(std::complex<float>)*samples_per_profile_);
          
          // Copy trajectory for profile
          //

          memcpy( trajectory->get_data_ptr() + 
                  frame*samples_per_profile_*profiles_per_frame_ + 
                  profile*samples_per_profile_,
                  
                  host_traj->get_data_ptr() + 
                  (first_profile + profile) * samples_per_profile_,
                  
                  sizeof(floatd2)*samples_per_profile_);
        }
      }
    }
    return GADGET_OK;
  }
  
  GadgetContainerMessage< hoNDArray< std::complex<float> > >*
  gpuRetroGatedSensePrepGadget::duplicate_profile( GadgetContainerMessage< hoNDArray< std::complex<float> > > *profile )
  {
    GadgetContainerMessage< hoNDArray< std::complex<float> > > *copy = 
      new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
    
    *copy->getObjectPtr() = *profile->getObjectPtr();
    
    return copy;
  }

  void gpuRetroGatedSensePrepGadget::reconfigure(unsigned int set, unsigned int slice)
  {    
    calculate_density_compensation_for_reconstruction(set, slice);
    
    cuSenseBuffer<float,2> *acc_buffer = (buffer_using_solver_) ? &acc_buffer_cg_[set*slices_+slice] : &acc_buffer_[set*slices_+slice];

    acc_buffer->setup( from_std_vector<size_t,2>(image_dimensions_recon_), image_dimensions_recon_os_, 
                       kernel_width_, num_coils_[set*slices_+slice],                        
                       num_buffer_frames_outer_, num_buffer_frames_inner_ );
    
    boost::shared_ptr< cuNDArray<float> > device_weights = calculate_density_compensation_for_buffer(set, slice);
    acc_buffer->set_dcw(device_weights);

    if( buffer_using_solver_ ){
      ((cuSenseBufferCg<float,2>*) acc_buffer)->set_dcw_for_rhs(calculate_density_compensation_for_rhs(set, slice));
      ((cuSenseBufferCg<float,2>*) acc_buffer)->preprocess(calculate_trajectory_for_rhs(0, set, slice).get());
    }
    
    reconfigure_[set*slices_+slice] = false;
  }

  GADGET_FACTORY_DECLARE(gpuRetroGatedSensePrepGadget)
}
