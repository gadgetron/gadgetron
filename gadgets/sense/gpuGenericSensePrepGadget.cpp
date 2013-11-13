#include "gpuGenericSensePrepGadget.h"
#include "Gadgetron.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "cuNonCartesianSenseOperator.h"
#include "SenseJob.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_utils.h"
#include "hoNDArray_utils.h"
#include "vector_td_operators.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "check_CUDA.h"
#include "hoNDArray_fileio.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace Gadgetron{

  gpuGenericSensePrepGadget::gpuGenericSensePrepGadget()
    : slices_(-1)
    , sets_(-1)
    , device_number_(-1)
    , samples_per_readout_(-1)
  {
    // Set some default values in case the config does not contain a specification
    //

    set_parameter(std::string("deviceno").c_str(), "0");
    set_parameter(std::string("rotations_per_reconstruction").c_str(), "0");
    set_parameter(std::string("propagate_csm_from_set").c_str(), "-1");
    set_parameter(std::string("buffer_length_in_rotations").c_str(), "0");
    set_parameter(std::string("buffer_using_solver").c_str(), "false");
    set_parameter(std::string("buffer_convolution_kernel_width").c_str(), "5.5");
    set_parameter(std::string("buffer_convolution_oversampling_factor").c_str(), "1.25");
    set_parameter(std::string("reconstruction_os_factor_x").c_str(), "1.0");
    set_parameter(std::string("reconstruction_os_factor_y").c_str(), "1.0");
  }
  
  gpuGenericSensePrepGadget::~gpuGenericSensePrepGadget() {}
  
  int gpuGenericSensePrepGadget::process_config(ACE_Message_Block* mb)
  {
    // Get configuration values from config file
    //

    device_number_ = get_int_value(std::string("deviceno").c_str());
    rotations_per_reconstruction_ = get_int_value(std::string("rotations_per_reconstruction").c_str());
    buffer_length_in_rotations_ = get_int_value(std::string("buffer_length_in_rotations").c_str());
    buffer_using_solver_ = get_bool_value(std::string("buffer_using_solver").c_str());
    output_timing_ = get_bool_value(std::string("output_timing").c_str());

    // Currently there are some restrictions on the allowed sliding window configurations
    //
    
    sliding_window_readouts_ = get_int_value(std::string("sliding_window_readouts").c_str());
    sliding_window_rotations_ = get_int_value(std::string("sliding_window_rotations").c_str());

    if( sliding_window_readouts_>0 && sliding_window_rotations_>0 ){
      GADGET_DEBUG1( "Error: Sliding window reconstruction is not yet supported for both readouts and frames simultaneously.\n" );
      return GADGET_FAIL;
    }

    if( sliding_window_readouts_>0 && rotations_per_reconstruction_>0 ){
      GADGET_DEBUG1( "Error: Sliding window reconstruction over readouts is not yet supported for multiframe reconstructions.\n" );
      return GADGET_FAIL;
    }
    
    if( sliding_window_rotations_ > 0 && sliding_window_rotations_ >= rotations_per_reconstruction_ ){
      GADGET_DEBUG1( "Error: Illegal sliding window configuration.\n" );
      return GADGET_FAIL;
    }

    // Setup and validate device configuration
    //

    int number_of_devices;
    if (cudaGetDeviceCount(&number_of_devices)!= cudaSuccess) {
      GADGET_DEBUG1( "Error: unable to query number of CUDA devices.\n" );
      return GADGET_FAIL;
    }

    if (number_of_devices == 0) {
      GADGET_DEBUG1( "Error: No available CUDA devices.\n" );
      return GADGET_FAIL;
    }

    if (device_number_ >= number_of_devices) {
      GADGET_DEBUG2("Adjusting device number from %d to %d\n", device_number_,  (device_number_%number_of_devices));
      device_number_ = (device_number_%number_of_devices);
    }

    if (cudaSetDevice(device_number_)!= cudaSuccess) {
      GADGET_DEBUG1( "Error: unable to set CUDA device.\n" );
      return GADGET_FAIL;
    }

    cudaDeviceProp deviceProp;
    if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
      GADGET_DEBUG1( "Error: unable to query device properties.\n" );
      return GADGET_FAIL;
    }
    
    unsigned int warp_size = deviceProp.warpSize;

    // It is possible to specify one set to use for csm propagation, and then propagate this to all sets
    //

    propagate_csm_from_set_ = get_int_value(std::string("propagate_csm_from_set").c_str());

    if( propagate_csm_from_set_ > 0 ){
      GADGET_DEBUG2("Currently, only set 0 can propagate coil sensitivity maps. Set %d was specified.\n", propagate_csm_from_set_ );
      return GADGET_FAIL;
    }

    if( propagate_csm_from_set_ >= 0 ){
      GADGET_DEBUG2("Propagating csm from set %d to all sets\n", propagate_csm_from_set_ );
    }

    // Convolution kernel width and oversampling ratio (for the buffer)
    //

    kernel_width_ = get_double_value(std::string("buffer_convolution_kernel_width").c_str());
    oversampling_factor_ = get_double_value(std::string("buffer_convolution_oversampling_factor").c_str());

    // Get the Ismrmrd header
    //

    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));
    
    if( cfg.get() == 0x0 ){
      GADGET_DEBUG1("Unable to parse Ismrmrd header\n");
      return GADGET_FAIL;
    }

    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();

    if (e_seq.size() != 1) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
      GADGET_DEBUG1("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }
    
    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    // Matrix sizes (as a multiple of the GPU's warp size)
    //
    
    image_dimensions_.push_back(e_space.matrixSize().x());
    image_dimensions_.push_back(e_space.matrixSize().y());

    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(e_space.matrixSize().x()*get_double_value(std::string("reconstruction_os_factor_x").c_str())))+warp_size-1)/warp_size)*warp_size);  

    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(e_space.matrixSize().y()*get_double_value(std::string("reconstruction_os_factor_y").c_str())))+warp_size-1)/warp_size)*warp_size);
    
    image_dimensions_recon_os_ = uintd2
      (((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
       ((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);
    
    // In case the warp_size constraint kicked in
    oversampling_factor_ = float(image_dimensions_recon_os_[0])/float(image_dimensions_recon_[0]); 
    
    GADGET_DEBUG2("matrix_size_x : %d, recon: %d, recon_os: %d\n", 
                  image_dimensions_[0], image_dimensions_recon_[0], image_dimensions_recon_os_[0]);

    GADGET_DEBUG2("matrix_size_y : %d, recon: %d, recon_os: %d\n", 
                  image_dimensions_[1], image_dimensions_recon_[1], image_dimensions_recon_os_[1]);
    
    fov_.push_back(r_space.fieldOfView_mm().x());
    fov_.push_back(r_space.fieldOfView_mm().y());
    fov_.push_back(r_space.fieldOfView_mm().z());

    slices_ = e_limits.slice().present() ? e_limits.slice().get().maximum() + 1 : 1;
    sets_ = e_limits.set().present() ? e_limits.set().get().maximum() + 1 : 1;
    
    // Allocate readout and trajectory queues
    // - one queue for the currently incoming frame
    // - one queue for the upcoming reconstruction

    frame_readout_queue_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);
    recon_readout_queue_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);
    frame_traj_queue_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);
    recon_traj_queue_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);
    image_headers_queue_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);
    
    size_t bsize = sizeof(GadgetContainerMessage< hoNDArray< std::complex<float> > >)*image_dimensions_[0]*10;
    
    for( unsigned int i=0; i<slices_*sets_; i++ ){
      frame_readout_queue_[i].high_water_mark(bsize);
      frame_readout_queue_[i].low_water_mark(bsize);
      frame_traj_queue_[i].high_water_mark(bsize);
      frame_traj_queue_[i].low_water_mark(bsize);
    }
    
    bsize *= (rotations_per_reconstruction_+1);
    
    for( unsigned int i=0; i<slices_*sets_; i++ ){
      recon_readout_queue_[i].high_water_mark(bsize);
      recon_readout_queue_[i].low_water_mark(bsize);
      recon_traj_queue_[i].high_water_mark(bsize);
      recon_traj_queue_[i].low_water_mark(bsize);
    }
    
    // Define various per slice/set variables
    //

    previous_readout_no_ = boost::shared_array<long>(new long[slices_*sets_]);
    acceleration_factor_ = boost::shared_array<long>(new long[slices_*sets_]);
    image_counter_ = boost::shared_array<long>(new long[slices_*sets_]);
    readout_counter_frame_= boost::shared_array<long>(new long[slices_*sets_]);
    readout_counter_global_= boost::shared_array<long>(new long[slices_*sets_]);
    readouts_per_frame_= boost::shared_array<long>(new long[slices_*sets_]);
    frames_per_rotation_= boost::shared_array<long>(new long[slices_*sets_]);
    buffer_frames_per_rotation_= boost::shared_array<long>(new long[slices_*sets_]);
    buffer_update_needed_ = boost::shared_array<bool>(new bool[slices_*sets_]);
    reconfigure_ = boost::shared_array<bool>(new bool[slices_*sets_]);
    num_coils_ = boost::shared_array<unsigned int>(new unsigned int[slices_*sets_]);
    
    for( unsigned int i=0; i<slices_*sets_; i++ ){

      previous_readout_no_[i] = -1;
      acceleration_factor_[i] = -1;
      image_counter_[i] = 0;
      readout_counter_frame_[i] = 0;
      readout_counter_global_[i] = 0;
      readouts_per_frame_[i] = get_int_value(std::string("readouts_per_frame").c_str());
      frames_per_rotation_[i] = get_int_value(std::string("frames_per_rotation").c_str());
      buffer_frames_per_rotation_[i] = get_int_value(std::string("buffer_frames_per_rotation").c_str());
      num_coils_[i] = 0;
      buffer_update_needed_[i] = true;
      reconfigure_[i] = true;

      // Assign some default values ("upper bound estimates") of the (possibly) unknown entities
      //
      
      if( readouts_per_frame_[i] == 0 ){
        readouts_per_frame_[i] = image_dimensions_[0];
      }
      
      if( frames_per_rotation_[i] == 0 ){
        frames_per_rotation_[i] = image_dimensions_[0]/readouts_per_frame_[i];
      }

      // Also remember to set the high/low water marks of the ISMRMRD image header queue
      //

      bsize = sizeof(GadgetContainerMessage<ISMRMRD::ImageHeader>)*100*
        std::max(1L, frames_per_rotation_[i]*rotations_per_reconstruction_);
    
      image_headers_queue_[i].high_water_mark(bsize);
      image_headers_queue_[i].low_water_mark(bsize);
    }

    // If need be the following limitation can be lifted, but it would be a little tedious... 
    //

    if( buffer_using_solver_ && rotations_per_reconstruction_ < 1 ) {
      GADGET_DEBUG1("Error: when buffering using a cg solver, 'rotations_per_reconstruction' must be specified (and strictly positive).");
    }

    if( buffer_using_solver_ && ( buffer_frames_per_rotation_[0] > 0 || buffer_length_in_rotations_ > 0 ) ){
      GADGET_DEBUG1("Error: when buffering using a cg solver, we currently do not support specification of 'buffer_frames_per_rotation' or 'buffer_length_in_rotations'. These values are instead automatically set to match the reconstruction settings.\n");
      return GADGET_FAIL;
    }
            
    position_ = boost::shared_array<float[3]>(new float[slices_*sets_][3]);
    read_dir_ = boost::shared_array<float[3]>(new float[slices_*sets_][3]);
    phase_dir_ = boost::shared_array<float[3]>(new float[slices_*sets_][3]);
    slice_dir_ = boost::shared_array<float[3]>(new float[slices_*sets_][3]);

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

    return GADGET_OK;
  }

  int gpuGenericSensePrepGadget::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,           // header
          GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2,   // data
          GadgetContainerMessage< hoNDArray<float> > *m3)                   // traj/dcw
  {
    // Noise should have been consumed by the noise adjust (if in the gadget chain)
    //
    
    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    if (is_noise) { 
      m1->release();
      return GADGET_OK;
    }

    // Setup timer if asked for
    //

    boost::shared_ptr<GPUTimer> process_timer;
    if( output_timing_ )
      process_timer = boost::shared_ptr<GPUTimer>( new GPUTimer("gpuGenericSensePrepGadget::process()") );

    // Some convienient utility variables
    //

    unsigned int set = m1->getObjectPtr()->idx.set;
    unsigned int slice = m1->getObjectPtr()->idx.slice;
    unsigned int readout = m1->getObjectPtr()->idx.kspace_encode_step_1;
    unsigned int idx = set*slices_+slice;

    // Get a pointer to the accumulation buffer. 
    //

    cuSenseBuffer<float,2> *acc_buffer = 
      (buffer_using_solver_) ? &acc_buffer_cg_[idx] : &acc_buffer_[idx];

    // Have the imaging plane changed?
    //

    if( !vec_equal(position_[idx], m1->getObjectPtr()->position) ||
        !vec_equal(read_dir_[idx], m1->getObjectPtr()->read_dir) || 
        !vec_equal(phase_dir_[idx], m1->getObjectPtr()->phase_dir) ||
        !vec_equal(slice_dir_[idx], m1->getObjectPtr()->slice_dir) ){
      
      // Yes indeed, clear the accumulation buffer and update structs
      //

      acc_buffer->clear();
      buffer_update_needed_[idx] = true;
      
      memcpy(position_[idx],m1->getObjectPtr()->position,3*sizeof(float));
      memcpy(read_dir_[idx],m1->getObjectPtr()->read_dir,3*sizeof(float));
      memcpy(phase_dir_[idx],m1->getObjectPtr()->phase_dir,3*sizeof(float));
      memcpy(slice_dir_[idx],m1->getObjectPtr()->slice_dir,3*sizeof(float));
    }
    
    // Only when the first readout arrives, do we know the #samples/readout
    //

    if( samples_per_readout_ == -1 )      
      samples_per_readout_ = m1->getObjectPtr()->number_of_samples;
    
    if( samples_per_readout_ != m1->getObjectPtr()->number_of_samples ){
      GADGET_DEBUG1("Unexpected change in the readout length\n");
      return GADGET_FAIL;
    }
    
    bool new_frame_detected = false;

    // Reconfigure at first pass
    // - or if the number of coil changes
    // - or if the reconfigure_ flag is set

    if( num_coils_[idx] != m1->getObjectPtr()->active_channels ){
      GADGET_DEBUG1("Reconfiguring (the number of coils changed)\n");
      num_coils_[idx] = m1->getObjectPtr()->active_channels;
      reconfigure(set, slice);
    }

    if( reconfigure_[idx] ){
      GADGET_DEBUG1("Reconfiguring (due to boolean indicator)\n");
      reconfigure(set, slice);
    }

    // Keep track of the incoming readout ids
    // - to determine the number of readouts per frame
    // - to determine the number of frames per rotation

    if (previous_readout_no_[idx] >= 0) {

      if ( readout > previous_readout_no_[idx]) { 
        // This is not the last readout in the frame.
        // Make an estimate of the acceleration factor
        //
	
        long tmp_accel = readout - previous_readout_no_[idx];

        if( acceleration_factor_[idx] != tmp_accel )
          GADGET_DEBUG2("Detected an acceleration factor of %d\n", tmp_accel);
	
        acceleration_factor_[idx] = tmp_accel;
      }
      else{ 

        // This is the first readout in a new frame
        //

        if( get_int_value(std::string("readouts_per_frame").c_str()) == 0 &&
            readout_counter_frame_[idx] > 0 &&
            readout_counter_frame_[idx] != readouts_per_frame_[idx] ){ 

          // A new acceleration factor is detected
          //

          GADGET_DEBUG1("Reconfiguring (acceleration factor changed)\n");

          new_frame_detected = true;
          readouts_per_frame_[idx] = readout_counter_frame_[idx];

          // Assume that #frames/rotation equals the acceleration factor
          // If not, or if we cannot deduce the acceleration factor from the difference
          // of two subsequent readout ids, then 'frames_per_rotation' have to be specified in the config...
          //
	    
          if( get_int_value(std::string("frames_per_rotation").c_str()) == 0 ) {
            frames_per_rotation_[idx] = acceleration_factor_[idx];
          }
          reconfigure(set, slice);
        }
      }
    }
    previous_readout_no_[idx] = readout;

    // Enqueue readout
    // - unless 'new_frame_detected', then the current readout does not belong to the current frame and we delay enqueing

    if( !new_frame_detected ) {
      
      // Memory handling is easier if we make copies for our internal queues
      frame_readout_queue_[idx].enqueue_tail(duplicate_array(m2));
      recon_readout_queue_[idx].enqueue_tail(duplicate_array(m2));
      frame_traj_queue_[idx].enqueue_tail(duplicate_array(m3));
      recon_traj_queue_[idx].enqueue_tail(duplicate_array(m3));
    }

    // If the readout is the last of a "true frame" (ignoring any sliding window readouts)
    // - then update the accumulation buffer

    bool is_last_readout_in_frame = (readout_counter_frame_[idx] == readouts_per_frame_[idx]-1);
    is_last_readout_in_frame |= new_frame_detected;

    cuNDArray<floatd2> traj;
    cuNDArray<float> dcw;
    
    if( is_last_readout_in_frame ){

      // Get ready to update the csm/regularization buffer
      //

      // Extract this frame's samples 
      //

      boost::shared_ptr< hoNDArray<float_complext> > host_samples = 
        extract_samples_from_queue( &frame_readout_queue_[idx], false, set, slice );
            
      cuNDArray<float_complext> samples( host_samples.get() );

      // Extract this frame's trajectory and dcw.
      //

      extract_trajectory_and_dcw_from_queue( &frame_traj_queue_[idx], false, set, slice, 
                                             samples_per_readout_*readouts_per_frame_[idx], 1,
                                             &traj, &dcw );

      // Scale dcw weights to the are of the oversampled recon matrix size
      float scale_factor = float(prod(image_dimensions_recon_os_))/asum(&dcw);
      dcw *= scale_factor;
      
      // Add this frame to the buffer
      //

      acc_buffer->set_dcw(boost::shared_ptr< cuNDArray<float> >(new cuNDArray<float>(&dcw)));
      buffer_update_needed_[idx] |= acc_buffer->add_frame_data( &samples, &traj );
    }

    // Are we ready to reconstruct (downstream)?
    //

    long readouts_per_reconstruction = readouts_per_frame_[idx];

    if( rotations_per_reconstruction_ > 0 )
      readouts_per_reconstruction *= (frames_per_rotation_[idx]*rotations_per_reconstruction_);
    
    bool is_last_readout_in_reconstruction = ( recon_readout_queue_[idx].message_count() == readouts_per_reconstruction );

    // Prepare the image header for this frame
    // - if this is indeed the last profile of a new frame
    // - or if we are about to reconstruct due to 'sliding_window_profiles_' > 0
    
    if( is_last_readout_in_frame || 
        (is_last_readout_in_reconstruction && image_headers_queue_[idx].message_count() == 0) ){
      
      GadgetContainerMessage<ISMRMRD::ImageHeader> *header = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
      ISMRMRD::AcquisitionHeader *base_head = m1->getObjectPtr();

      {
        // Initialize header to all zeroes (there is a few fields we do not set yet)
        ISMRMRD::ImageHeader tmp = {0};
        *(header->getObjectPtr()) = tmp;
      }

      header->getObjectPtr()->version = base_head->version;

      header->getObjectPtr()->matrix_size[0] = image_dimensions_recon_[0];
      header->getObjectPtr()->matrix_size[1] = image_dimensions_recon_[1];
      header->getObjectPtr()->matrix_size[2] = std::max(1L,frames_per_rotation_[idx]*rotations_per_reconstruction_);

      header->getObjectPtr()->field_of_view[0] = fov_[0];
      header->getObjectPtr()->field_of_view[1] = fov_[1];
      header->getObjectPtr()->field_of_view[2] = fov_[2];

      header->getObjectPtr()->channels = num_coils_[idx];
      header->getObjectPtr()->slice = base_head->idx.slice;
      header->getObjectPtr()->set = base_head->idx.set;

      header->getObjectPtr()->acquisition_time_stamp = base_head->acquisition_time_stamp;
      memcpy(header->getObjectPtr()->physiology_time_stamp, base_head->physiology_time_stamp, sizeof(uint32_t)*ISMRMRD_PHYS_STAMPS);

      memcpy(header->getObjectPtr()->position, base_head->position, sizeof(float)*3);
      memcpy(header->getObjectPtr()->read_dir, base_head->read_dir, sizeof(float)*3);
      memcpy(header->getObjectPtr()->phase_dir, base_head->phase_dir, sizeof(float)*3);
      memcpy(header->getObjectPtr()->slice_dir, base_head->slice_dir, sizeof(float)*3);
      memcpy(header->getObjectPtr()->patient_table_position, base_head->patient_table_position, sizeof(float)*3);

      header->getObjectPtr()->image_data_type = ISMRMRD::DATA_COMPLEX_FLOAT;
      header->getObjectPtr()->image_index = image_counter_[idx]++; 
      header->getObjectPtr()->image_series_index = idx;

      image_headers_queue_[idx].enqueue_tail(header);
    }
    
    // If it is time to reconstruct (downstream) then prepare the Sense job
    // 

    if( is_last_readout_in_reconstruction ){
      
      // Update csm and regularization images if the buffer has changed (completed a cycle) 
      // - and at the first pass

      if( buffer_update_needed_[idx] || 
          csm_host_[idx].get_number_of_elements() == 0 || 
          reg_host_[idx].get_number_of_elements() == 0 ){

        // Get the accumulated coil images
        //
        
        boost::shared_ptr< cuNDArray<float_complext> > csm_data = acc_buffer->get_accumulated_coil_images();
        
        // Estimate CSM
        //
        
        if( propagate_csm_from_set_ < 0 || propagate_csm_from_set_ == set ){	  	  
          csm_ = estimate_b1_map<float,2>( csm_data.get() );
        }
        else{
          GADGET_DEBUG2("Set %d is reusing the csm from set %d\n", set, propagate_csm_from_set_);
          if( csm_.get() == 0x0 ){
            GADGET_DEBUG1("Error: csm has not been computed, cannot propagate\n");
            return GADGET_FAIL;
          }	  
        }

        acc_buffer->set_csm(csm_);
        csm_host_[idx] = *(csm_->to_host());
	
        // Compute regularization image
        //

        boost::shared_ptr< cuNDArray<float_complext> > reg_image;
        std::vector<unsigned int> dims;
    	
        if( buffer_using_solver_ ){

          //GPUTimer timer("\n\n AVOIDABLE PREPROCESSING. HOW EXPENSIVE?\n\n");

          extract_trajectory_and_dcw_from_queue( &recon_traj_queue_[idx], true, set, slice, 
                                                 samples_per_readout_*readouts_per_frame_[idx],
                                                 std::max(1L, frames_per_rotation_[idx]*rotations_per_reconstruction_),
                                                 &traj, &dcw );

          // Scale dcw weights to the are of the oversampled recon matrix size
          float scale_factor = float(prod(image_dimensions_recon_os_))/asum(&dcw);
          dcw *= scale_factor;

          dims = *traj.get_dimensions();

          std::vector<unsigned int> tmp_dims;
          tmp_dims.push_back(dims[0]*dims[1]);
          tmp_dims.push_back(1);
	  
          traj.reshape(&tmp_dims);
          dcw.reshape(&tmp_dims);
	  
          ((cuSenseBufferCg<float,2>*)acc_buffer)->preprocess(&traj);
          ((cuSenseBufferCg<float,2>*)acc_buffer)->set_dcw_for_rhs(boost::shared_ptr< cuNDArray<float> >(new cuNDArray<float>(&dcw)));
        }

        reg_image = acc_buffer->get_combined_coil_image();	
        reg_host_[idx] = *(reg_image->to_host());

        if( buffer_using_solver_ ){
          traj.reshape(&dims);
          dcw.reshape(&dims);
        }

        buffer_update_needed_[idx] = false;
      }

      // Prepare data array for the downstream reconstruction
      //
      
      boost::shared_ptr< hoNDArray<float_complext> > samples_host = 
        extract_samples_from_queue( &recon_readout_queue_[idx], true, set, slice );
      
      // Preapre the trajectory and dcw arrays.
      // They have already been computed above 
      // - if 'rotations_per_reconstruction_' is 0
      // - if 'buffer_using_solver_' is true
      
      if( !(/*rotations_per_reconstruction_ == 0 ||*/ buffer_using_solver_) ){
      	extract_trajectory_and_dcw_from_queue( &recon_traj_queue_[idx], true, set, slice, 
                                               samples_per_readout_*readouts_per_frame_[idx],
                                               std::max(1L, frames_per_rotation_[idx]*rotations_per_reconstruction_),
                                               &traj, &dcw );
      }

      // Set up the Sense job
      //

      GadgetContainerMessage< SenseJob > *sj = new GadgetContainerMessage<SenseJob>();
      	
      sj->getObjectPtr()->dat_host_ = samples_host;      
      sj->getObjectPtr()->tra_host_ = traj.to_host();
      sj->getObjectPtr()->dcw_host_ = dcw.to_host();
      sj->getObjectPtr()->csm_host_ = boost::shared_ptr< hoNDArray<float_complext> >( new hoNDArray<float_complext>(csm_host_[idx]));
      sj->getObjectPtr()->reg_host_ = boost::shared_ptr< hoNDArray<float_complext> >( new hoNDArray<float_complext>(reg_host_[idx]));
      
      // Pull the image headers out of the queue
      //

      long frames_per_reconstruction = 
        std::max( 1L, frames_per_rotation_[idx]*rotations_per_reconstruction_ );
      
      if( image_headers_queue_[idx].message_count() != frames_per_reconstruction ){
        sj->release();
        GADGET_DEBUG2("Unexpected size of image header queue: %d, %d\n", 
                      image_headers_queue_[idx].message_count(), frames_per_reconstruction);
        return GADGET_FAIL;
      }
      
      sj->getObjectPtr()->image_headers_ =
        boost::shared_array<ISMRMRD::ImageHeader>( new ISMRMRD::ImageHeader[frames_per_reconstruction] );
      
      for( unsigned int i=0; i<frames_per_reconstruction; i++ ){	

        ACE_Message_Block *mbq;

        if( image_headers_queue_[idx].dequeue_head(mbq) < 0 ) {
          sj->release();
          GADGET_DEBUG1("Image header dequeue failed\n");
          return GADGET_FAIL;
        }
	
        GadgetContainerMessage<ISMRMRD::ImageHeader> *m = AsContainerMessage<ISMRMRD::ImageHeader>(mbq);
        sj->getObjectPtr()->image_headers_[i] = *m->getObjectPtr();

        // In sliding window mode the header might need to go back at the end of the queue for reuse
        // 
	
        if( i >= frames_per_reconstruction-sliding_window_rotations_*frames_per_rotation_[idx] ){
          image_headers_queue_[idx].enqueue_tail(m);
        }
        else {
          m->release();
        }
      }
      
      // The Sense Job needs an image header as well. 
      // Let us just copy the initial one...

      GadgetContainerMessage<ISMRMRD::ImageHeader> *m4 = new GadgetContainerMessage<ISMRMRD::ImageHeader>;

      *m4->getObjectPtr() = sj->getObjectPtr()->image_headers_[0];
      m4->cont(sj);

      // Pass the Sense job downstream
      //
      
      if (this->next()->putq(m4) < 0) {
        GADGET_DEBUG1("Failed to put job on queue.\n");
        m4->release();
        return GADGET_FAIL;
      }
    }
    
    if( is_last_readout_in_frame )
      readout_counter_frame_[idx] = 0;
    else{
      readout_counter_frame_[idx]++;
    }

    if( new_frame_detected ){

      // The incoming profile was actually the first readout of the next frame, enqueue.
      //

      frame_readout_queue_[idx].enqueue_tail(duplicate_array(m2));
      recon_readout_queue_[idx].enqueue_tail(duplicate_array(m2)); 
      frame_traj_queue_[idx].enqueue_tail(duplicate_array(m3));
      recon_traj_queue_[idx].enqueue_tail(duplicate_array(m3)); 

      readout_counter_frame_[idx]++;
    }

    readout_counter_global_[idx]++;

    if( output_timing_ )
      process_timer.reset();
    
    m1->release(); // this is safe, the internal queues hold copies
    return GADGET_OK;
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > 
  gpuGenericSensePrepGadget::extract_samples_from_queue ( ACE_Message_Queue<ACE_MT_SYNCH> *queue, 
                                                          bool sliding_window, unsigned int set, unsigned int slice )
  {    
    unsigned int readouts_buffered = queue->message_count();
    
    std::vector<unsigned int> dims;
    dims.push_back(samples_per_readout_*readouts_buffered);
    dims.push_back(num_coils_[set*slices_+slice]);
    
    boost::shared_ptr< hoNDArray<float_complext> > host_samples(new hoNDArray<float_complext>(&dims));
    
    for (unsigned int p=0; p<readouts_buffered; p++) {
      
      ACE_Message_Block* mbq;
      if (queue->dequeue_head(mbq) < 0) {
        GADGET_DEBUG1("Message dequeue failed\n");
        throw std::runtime_error("gpuGenericSensePrepGadget::extract_samples_from_queue: dequeing failed");	
      }
      
      GadgetContainerMessage< hoNDArray< std::complex<float> > > *daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq);
	
      if (!daq) {
        GADGET_DEBUG1("Unable to interpret data on message queue\n");
        throw std::runtime_error("gpuGenericSensePrepGadget::extract_samples_from_queue: failed to interpret data");	
      }
	
      for (unsigned int c = 0; c < num_coils_[set*slices_+slice]; c++) {
	
        float_complext *data_ptr = host_samples->get_data_ptr();
        data_ptr += c*samples_per_readout_*readouts_buffered+p*samples_per_readout_;
	    
        std::complex<float> *r_ptr = daq->getObjectPtr()->get_data_ptr();
        r_ptr += c*daq->getObjectPtr()->get_size(0);
	  
        memcpy(data_ptr, r_ptr, samples_per_readout_*sizeof(float_complext));
      }

      // In sliding window mode the readout might need to go back at the end of the queue
      // 
      
      long readouts_in_sliding_window = sliding_window_readouts_ + 
        readouts_per_frame_[set*slices_+slice]*frames_per_rotation_[set*slices_+slice]*sliding_window_rotations_;

      if( sliding_window && p >= (readouts_buffered-readouts_in_sliding_window) )
        queue->enqueue_tail(mbq);
      else
        mbq->release();
    } 
    
    return host_samples;
  }
  
  boost::shared_ptr< hoNDArray<float> > 
  gpuGenericSensePrepGadget::extract_trajectory_from_queue ( ACE_Message_Queue<ACE_MT_SYNCH> *queue, 
                                                             bool sliding_window, unsigned int set, unsigned int slice )
  {    
    if(!queue) {
      GADGET_DEBUG1("Illegal queue pointer, cannot extract trajectory\n");
      throw std::runtime_error("gpuGenericSensePrepGadget::extract_trajectory_from_queue: illegal queue pointer");	
    }

    if(queue->message_count()==0) {
      GADGET_DEBUG1("Empty queue, cannot extract trajectory\n");
      throw std::runtime_error("gpuGenericSensePrepGadget::extract_trajectory_from_queue: empty queue");	
    }

    if(samples_per_readout_ < 1) {
      GADGET_DEBUG2("Empty queue (%d), cannot extract trajectory\n", samples_per_readout_);
      throw std::runtime_error("gpuGenericSensePrepGadget::extract_trajectory_from_queue: empty queue");	
    }
    
    unsigned int readouts_buffered = queue->message_count();
    
    std::vector<unsigned int> dims;
    dims.push_back(3);
    dims.push_back(samples_per_readout_);
    dims.push_back(readouts_buffered);
    
    boost::shared_ptr< hoNDArray<float> > host_samples(new hoNDArray<float>(&dims));
    
    for (unsigned int p=0; p<readouts_buffered; p++) {      
      ACE_Message_Block* mbq;
      if (queue->dequeue_head(mbq) < 0) {
        GADGET_DEBUG1("Message dequeue failed\n");
        throw std::runtime_error("gpuGenericSensePrepGadget::extract_trajectory_from_queue: dequeing failed");	
      }
      
      GadgetContainerMessage< hoNDArray<float> > *daq = AsContainerMessage<hoNDArray<float> >(mbq);
	
      if (!daq) {
        GADGET_DEBUG1("Unable to interpret data on message queue\n");
        throw std::runtime_error("gpuGenericSensePrepGadget::extract_trajectory_from_queue: failed to interpret data");	
      }

      float *data_ptr = host_samples->get_data_ptr();
      data_ptr += 3*samples_per_readout_*p;
      
      float *r_ptr = daq->getObjectPtr()->get_data_ptr();
      
      memcpy(data_ptr, r_ptr, 3*samples_per_readout_*sizeof(float));
      
      // In sliding window mode the readout might need to go back at the end of the queue
      // 
      
      long readouts_in_sliding_window = sliding_window_readouts_ + 
        readouts_per_frame_[set*slices_+slice]*frames_per_rotation_[set*slices_+slice]*sliding_window_rotations_;

      if( sliding_window && p >= (readouts_buffered-readouts_in_sliding_window) )
        queue->enqueue_tail(mbq);
      else
        mbq->release();
    } 
    
    return host_samples;
  }
  
  void gpuGenericSensePrepGadget::extract_trajectory_and_dcw_from_queue
  ( ACE_Message_Queue<ACE_MT_SYNCH> *queue, bool sliding_window, unsigned int set, unsigned int slice, 
    unsigned int samples_per_frame, unsigned int num_frames,
    cuNDArray<floatd2> *traj, cuNDArray<float> *dcw )
  {
    // Extract trajectory and dcw.
    // They are stored as a float array of dimensions: 3 x #samples_per_readout x #readouts.
    // We need
    // - a floatd2 trajectory array 
    // - a float dcw array 
    //
    
    boost::shared_ptr< hoNDArray<float> > host_traj_dcw =
      extract_trajectory_from_queue( queue, sliding_window, set, slice );
    
    std::vector<unsigned int> order;
    order.push_back(1); order.push_back(2); order.push_back(0);
    
    boost::shared_ptr< hoNDArray<float> > host_traj_dcw_shifted =
      permute( host_traj_dcw.get(), &order );
    
    std::vector<unsigned int> dims_1d;
    dims_1d.push_back(host_traj_dcw_shifted->get_size(0)*host_traj_dcw_shifted->get_size(1));
    
    {
      hoNDArray<float> tmp(&dims_1d, host_traj_dcw_shifted->get_data_ptr()+2*dims_1d[0]);
      *dcw = tmp;
    }
    
    std::vector<unsigned int> dims_2d = dims_1d;
    dims_2d.push_back(2);
    
    order.clear();
    order.push_back(1); order.push_back(0);

    hoNDArray<float> tmp(&dims_2d, host_traj_dcw_shifted->get_data_ptr());
    cuNDArray<float> __traj(&tmp);
    boost::shared_ptr< cuNDArray<float> > _traj = permute( &__traj, &order );
    
    cuNDArray<floatd2> tmp2(&dims_1d, (floatd2*)_traj->get_data_ptr());
    
    *traj = tmp2;
    
    unsigned int idx = set*slices_+slice;
    dims_2d.clear();

    dims_2d.push_back(samples_per_frame);
    dims_2d.push_back(num_frames);

    dcw->reshape(&dims_2d);
    traj->reshape(&dims_2d);
  }

  template<class T> GadgetContainerMessage< hoNDArray<T> >*
  gpuGenericSensePrepGadget::duplicate_array( GadgetContainerMessage< hoNDArray<T> > *array )
  {
    GadgetContainerMessage< hoNDArray<T> > *copy = new GadgetContainerMessage< hoNDArray<T> >();   
    *(copy->getObjectPtr()) = *(array->getObjectPtr());
    return copy;
  }

  void gpuGenericSensePrepGadget::reconfigure(unsigned int set, unsigned int slice)
  {    
    unsigned int idx = set*slices_+slice;
    
    GADGET_DEBUG2("\nReconfiguring:\n#readouts/frame:%d\n#frames/rotation: %d\n#rotations/reconstruction:%d\n", 
                  readouts_per_frame_[idx], frames_per_rotation_[idx], rotations_per_reconstruction_);
    
    buffer_frames_per_rotation_[idx] = get_int_value(std::string("buffer_frames_per_rotation").c_str());
    
    if( buffer_frames_per_rotation_[idx] == 0 ){
      buffer_frames_per_rotation_[idx] = frames_per_rotation_[idx];
    }
    
    if( get_int_value(std::string("buffer_length_in_rotations").c_str()) == 0 ){
      buffer_length_in_rotations_ = std::max(1L, rotations_per_reconstruction_);
    }

    cuSenseBuffer<float,2> *acc_buffer = 
      (buffer_using_solver_) ? &acc_buffer_cg_[idx] : &acc_buffer_[idx];
    
    if( buffer_frames_per_rotation_[idx] == 1 ){ // Is this general enough to detect golden ratio type trajectories?

      acc_buffer->setup( from_std_vector<unsigned int,2>(image_dimensions_recon_), image_dimensions_recon_os_, 
                         kernel_width_, num_coils_[idx], 1, buffer_length_in_rotations_ );
    }else{
      acc_buffer->setup( from_std_vector<unsigned int,2>(image_dimensions_recon_), image_dimensions_recon_os_, 
                         kernel_width_, num_coils_[idx], buffer_length_in_rotations_, buffer_frames_per_rotation_[idx] );
    }
    reconfigure_[idx] = false;
  }

  GADGET_FACTORY_DECLARE(gpuGenericSensePrepGadget)
}
