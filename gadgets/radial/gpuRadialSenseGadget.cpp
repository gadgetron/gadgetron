#include "gpuRadialSenseGadget.h"
#include "Gadgetron.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "SenseJob.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_utils.h"
#include "vector_td_operators.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "check_CUDA.h"
#include "radial_utilities.h"
#include "hoNDArray_fileio.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <boost/shared_ptr.hpp>

namespace Gadgetron{

  gpuRadialSenseGadget::gpuRadialSenseGadget()
    : slices_(-1)
    , sets_(-1)
    , image_counter_(0)
    , image_series_(0)
    , device_number_(-1)
    , mode_(-1)
    , samples_per_profile_(-1)
    , previous_profile_(0x0)
    , profiles_counter_frame_(0x0)
    , profiles_counter_global_(0x0)
    , kernel_width_(-1.0f)
    , oversampling_factor_(-1.0f)
    , buffer_update_needed_(false)
    , reconfigure_(true)
  {
    // Set some default values in case the config does not contain a specification
    //

    set_parameter(std::string("mode").c_str(), "0");
    set_parameter(std::string("deviceno").c_str(), "0");
    set_parameter(std::string("buffer_using_solver").c_str(), "false");
    set_parameter(std::string("buffer_convolution_kernel_width").c_str(), "5.5");
    set_parameter(std::string("buffer_convolution_oversampling_factor").c_str(), "1.25");
    set_parameter(std::string("rotations_per_reconstruction").c_str(), "0");
    set_parameter(std::string("reconstruction_os_factor_x").c_str(), "1.0");
    set_parameter(std::string("reconstruction_os_factor_y").c_str(), "1.0");
  }
  
  gpuRadialSenseGadget::~gpuRadialSenseGadget()
  {
    if( previous_profile_ ) delete[] previous_profile_;
    if( profiles_counter_frame_ ) delete[] profiles_counter_frame_;
    if( profiles_counter_global_ ) delete[] profiles_counter_global_;
  }
  
  int gpuRadialSenseGadget::process_config(ACE_Message_Block* mb)
  {
    //GADGET_DEBUG1("gpuRadialSenseGadget::process_config\n");

    mode_ = get_int_value(std::string("mode").c_str());
    device_number_ = get_int_value(std::string("deviceno").c_str());
    profiles_per_frame_ = get_int_value(std::string("profiles_per_frame").c_str());
    frames_per_rotation_ = get_int_value(std::string("frames_per_rotation").c_str());
    rotations_per_reconstruction_ = get_int_value(std::string("rotations_per_reconstruction").c_str());
    buffer_length_in_rotations_ = get_int_value(std::string("buffer_length_in_rotations").c_str());
    buffer_using_solver_ = get_bool_value(std::string("buffer_using_solver").c_str());

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

    // Convolution kernel width and oversampling ratio (for the buffer)
    kernel_width_ = get_double_value(std::string("buffer_convolution_kernel_width").c_str());
    oversampling_factor_ = get_double_value(std::string("buffer_convolution_oversampling_factor").c_str());

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
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    // Matrix sizes (as a multiple of the GPU's warp size)

    image_dimensions_.push_back(((e_space.matrixSize().x()+warp_size-1)/warp_size)*warp_size);
    image_dimensions_.push_back(((e_space.matrixSize().y()+warp_size-1)/warp_size)*warp_size);
    
    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(e_space.matrixSize().x()*get_double_value(std::string("reconstruction_os_factor_x").c_str())))+warp_size-1)/warp_size)*warp_size);  
    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(e_space.matrixSize().y()*get_double_value(std::string("reconstruction_os_factor_y").c_str())))+warp_size-1)/warp_size)*warp_size);
    
    uintd2 matrix_size = uintd2(image_dimensions_recon_[0],image_dimensions_recon_[1]);
    image_dimensions_recon_os_ = uintd2
      (((static_cast<unsigned int>(std::ceil(matrix_size.vec[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
       ((static_cast<unsigned int>(std::ceil(matrix_size.vec[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);
    
    oversampling_factor_ = float(image_dimensions_recon_os_[0])/float(matrix_size[0]); // in case the warp_size constraint kicked in
    
    GADGET_DEBUG2("matrix_size_x : %d, recon: %d, recon_os: %d\n", image_dimensions_[0], image_dimensions_recon_[0], image_dimensions_recon_os_[0]);
    GADGET_DEBUG2("matrix_size_y : %d, recon: %d, recon_os: %d\n", image_dimensions_[1], image_dimensions_recon_[1], image_dimensions_recon_os_[1]);
    
    slices_ = e_limits.slice().present() ? e_limits.slice().get().maximum() + 1 : 1;
    sets_ = e_limits.set().present() ? e_limits.set().get().maximum() + 1 : 1;
    
    // Allocate profile queues
    // - one queue for the currently incoming frame
    // - one queue for the next reconstruction

    frame_queue_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);
    recon_queue_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);

    size_t bsize = (sizeof(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>)+
		    sizeof(GadgetContainerMessage< hoNDArray< std::complex<float> > >))*image_dimensions_[0]*10;
  
    for( unsigned int i=0; i<slices_*sets_; i++ ){
      frame_queue_[i].high_water_mark(bsize);
      frame_queue_[i].low_water_mark(bsize);
    }
    
    bsize *= (rotations_per_reconstruction_+1);
    
    for( unsigned int i=0; i<slices_*sets_; i++ ){
      recon_queue_[i].high_water_mark(bsize);
      recon_queue_[i].low_water_mark(bsize);
      }
    
    previous_profile_ = new long[slices_*sets_];
    profiles_counter_frame_= new long[slices_*sets_];
    profiles_counter_global_= new long[slices_*sets_];

    if( !previous_profile_ || !profiles_counter_frame_ || !profiles_counter_global_ ){
      GADGET_DEBUG1("Failed to allocate host memory\n");
      return GADGET_FAIL;
    }
    
    for( unsigned int i=0; i<slices_*sets_; i++ ){
      previous_profile_[i] = -1;
      profiles_counter_frame_[i] = 0;
      profiles_counter_global_[i] = 0;
    }
    
    // Assign some default values ("upper bound estimates") of the yet unknown entities
    
    if( profiles_per_frame_ == 0 ){
      profiles_per_frame_ = image_dimensions_[0];
    }
    
    if( frames_per_rotation_ == 0 ){
      if( mode_ == 2 ) // golden angle
	frames_per_rotation_ = 1;
      else
	frames_per_rotation_ = image_dimensions_[0]/profiles_per_frame_;
    }

    // Allocate accumulation buffer
    if( buffer_using_solver_ )
      acc_buffer_ = boost::shared_ptr< cuSenseBufferCg<float,2> >(new cuSenseBufferCg<float,2>());
    else
      acc_buffer_ = boost::shared_ptr< cuSenseBuffer<float,2> >(new cuSenseBuffer<float,2>());
    
    return GADGET_OK;
  }

  int gpuRadialSenseGadget::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    // Noise should have been consumed by the noise adjust, but just in case...
    //
    
    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    if (is_noise) { 
      m1->release();
      return GADGET_OK;
    }

    //GADGET_DEBUG1("gpuRadialSenseGadget::process\n");
    //GPUTimer timer("gpuRadialGadget::process");

    // Have the imaging plane changed?
    //

    if ( !read_dir_equal(m1->getObjectPtr()->read_dir) || 
	 !phase_dir_equal(m1->getObjectPtr()->phase_dir) ||
	 !slice_dir_equal(m1->getObjectPtr()->slice_dir) ||
	 !position_equal(m1->getObjectPtr()->position)) {

      // If so then clear the accumulation buffer
      acc_buffer_->clear();
      
      memcpy(position_,m1->getObjectPtr()->position,3*sizeof(float));
      memcpy(read_dir_,m1->getObjectPtr()->read_dir,3*sizeof(float));
      memcpy(phase_dir_,m1->getObjectPtr()->phase_dir,3*sizeof(float));
      memcpy(slice_dir_,m1->getObjectPtr()->slice_dir,3*sizeof(float));
    }
    
    // Only when the first profile arrives, do we know the #samples/profile
    //

    if( samples_per_profile_ == -1 )      
      samples_per_profile_ = m1->getObjectPtr()->number_of_samples;
    
    if( samples_per_profile_ != m1->getObjectPtr()->number_of_samples ){
      GADGET_DEBUG1("Unexpected change in the incoming profiles' lengths\n");
      return GADGET_FAIL;
    }
    
    unsigned int profile = m1->getObjectPtr()->idx.kspace_encode_step_1;
    unsigned int slice = m1->getObjectPtr()->idx.slice;
    unsigned int set = m1->getObjectPtr()->idx.set;
    unsigned int num_coils = m1->getObjectPtr()->active_channels;
    bool new_frame_detected = false;

    // Reconfigure at first pass (and if the reconfigure_ flag is set later on)
    // 

    if( reconfigure_ )
      reconfigure(num_coils);

    // Keep track of the incoming profile ids (mode dependent)
    // - to determine the number of profiles per frame
    // - to determine the number of frames per rotation
    //

    if (previous_profile_[set*slices_+slice] >= 0) {

      if ( profile > previous_profile_[set*slices_+slice]) { // this is not the last profile in the frame
	if( mode_ == 0 && get_int_value(std::string("frames_per_rotation").c_str()) == 0 ){
	  unsigned int acceleration_factor = profile - previous_profile_[set*slices_+slice];
	  if( acceleration_factor != frames_per_rotation_ ){
	    frames_per_rotation_ = acceleration_factor;
	    reconfigure(num_coils);
	  }
	}
      }
      else{ // This is the first profile in a new frame
	if( get_int_value(std::string("profiles_per_frame").c_str()) == 0 && // make sure the user did not specify a desired value for this variable
	    profiles_counter_frame_[set*slices_+slice] > 0 &&
	    profiles_counter_frame_[set*slices_+slice] != profiles_per_frame_ ){ // a new acceleration factor is detected
	  new_frame_detected = true;
	  profiles_per_frame_ = profiles_counter_frame_[set*slices_+slice];
	  if( mode_ == 1 && get_int_value(std::string("frames_per_rotation").c_str()) == 0 )
	    frames_per_rotation_ = image_dimensions_[0]/profiles_per_frame_;
	  reconfigure(num_coils);
	}
      }
    }
    previous_profile_[set*slices_+slice] = profile;

    // Are we ready to add a frame to the buffer? Or reconstruct downstream?
    //

    long profiles_per_reconstruction = profiles_per_frame_*frames_per_rotation_*rotations_per_reconstruction_;

    bool is_last_profile_in_frame = (profiles_counter_frame_[set*slices_+slice] == profiles_per_frame_-1);
    is_last_profile_in_frame |= new_frame_detected;
       
    bool is_last_rotation = (rotations_per_reconstruction_ == 0);
    if( profiles_per_reconstruction > 0 )
      is_last_rotation |= ((profiles_counter_global_[set*slices_+slice]%profiles_per_reconstruction) == profiles_per_reconstruction-1);

    // Enqueue profile
    // - if 'new_frame_detected' the current profile does not belong to the current frame and we delay enqueing

    if( !new_frame_detected ) {
      
      // The memory handling is easier if we make copies for our internal queues
      frame_queue_[set*slices_+slice].enqueue_tail(duplicate_profile(m1));
      recon_queue_[set*slices_+slice].enqueue_tail(duplicate_profile(m1));
    }

    // If the profile is the last of the frame then update the accumulation buffer
    //

    if( is_last_profile_in_frame ){

      //GPUTimer timer("gpuRadialGadget::UPDATING BUFFER FRAME");

      std::vector<unsigned int> dims;
      dims.push_back(samples_per_profile_*profiles_per_frame_);
      dims.push_back(num_coils);

      boost::shared_ptr< hoNDArray<float_complext> > host_samples = 
	extract_samples_from_queue( &frame_queue_[set*slices_+slice], dims );

      if( host_samples.get() == 0x0 ){
	GADGET_DEBUG1("Failed to extract frame data from queue\n");
	return GADGET_FAIL;
      }
      
      cuNDArray<float_complext> samples( host_samples.get() );
      
      long profile_offset = profiles_counter_global_[set*slices_+slice] - ((new_frame_detected) ? 1 : 0);
      boost::shared_ptr< cuNDArray<floatd2> > traj = calculate_trajectory_for_frame(profile_offset);

      buffer_update_needed_ |= acc_buffer_->add_frame_data( &samples, traj.get() );
    }
    
    // If the profile is the last before a reconstruction then prepare the Sense job
    //
    
    if( is_last_profile_in_frame && is_last_rotation ){

      // Update csm and regularization images if the buffer has changed (completed a cycle) 
      // - and at the first pass

      if( buffer_update_needed_ || csm_host_.get() == 0x0 || reg_host_.get() == 0x0 ){

	//GPUTimer timer("gpuRadialGadget::UPDATING BUFFER CSM/REG");
      
	// Get the accumulated coil images
	//

	boost::shared_ptr< cuNDArray<float_complext> > csm_data = acc_buffer_->get_accumulated_coil_images();

	if( !csm_data.get() ){
	  GADGET_DEBUG1("Error during accumulation buffer computation\n");
	  return GADGET_FAIL;
	}            
	
	// Estimate CSM
	//

	boost::shared_ptr< cuNDArray<float_complext> > csm = estimate_b1_map<float,2>( csm_data.get() );

	if( !csm.get() ){
	  GADGET_DEBUG1("Error during coil estimation\n");
	  return GADGET_FAIL;
	}            

	csm_host_ = csm->to_host();
	
	// Compute regularization image
	//

	boost::shared_ptr< cuNDArray<float_complext> > reg_image;	
	acc_buffer_->set_csm(csm);
	
	if( buffer_using_solver_ && mode_ == 2 ){
	  cuSenseBufferCg<float,2> *acc = (cuSenseBufferCg<float,2>*) acc_buffer_.get();
	  acc->preprocess(calculate_trajectory_for_rhs(profiles_counter_global_[set*slices_+slice] - ((new_frame_detected) ? 1 : 0)).get());
	}

	reg_image = acc_buffer_->get_combined_coil_image();
	
	if( !reg_image.get() ){
	  GADGET_DEBUG1("Error computing regularization image\n");
	  return GADGET_FAIL;
	}            
	
	reg_host_ = reg_image->to_host();
		
	/*
	static int counter = 0;
	char filename[256];
	sprintf((char*)filename, "reg_%d.cplx", counter);
	write_nd_array<float_complext>( reg_host_.get(), filename );
	counter++; */

	buffer_update_needed_ = false;
      }

      // Prepare data array of the profiles for the downstream reconstruction
      //
      
      unsigned int num_frames = std::max(1L,frames_per_rotation_*rotations_per_reconstruction_);

      std::vector<unsigned int> dims;
      dims.push_back(samples_per_profile_*profiles_per_frame_*num_frames);
      dims.push_back(num_coils);

      boost::shared_ptr< hoNDArray<float_complext> > samples_host = 
	extract_samples_from_queue( &recon_queue_[set*slices_+slice], dims );
      
      if( samples_host.get() == 0x0 ){
	GADGET_DEBUG1("Failed to extract frame data from queue\n");
	return GADGET_FAIL;
      }
           
      // The trajectory needs to be updated on the fly:
      // - for golden ratio based acquisitions
      // - when we are reconstructing frame-by-frame
      
      if( mode_ == 2 || rotations_per_reconstruction_ == 0 ){
	calculate_trajectory_for_reconstruction( profiles_counter_global_[set*slices_+slice] - ((new_frame_detected) ? 1 : 0) );
      }
      
      // Set up Sense job
      //

      GadgetContainerMessage< SenseJob >* m4 = new GadgetContainerMessage< SenseJob >();
	
      m4->getObjectPtr()->dat_host_ = samples_host;
      m4->getObjectPtr()->tra_host_ = boost::shared_ptr< hoNDArray<floatd2> >(new hoNDArray<floatd2>(*host_traj_recon_));
      m4->getObjectPtr()->dcw_host_ = boost::shared_ptr< hoNDArray<float> >(new hoNDArray<float>(*host_weights_recon_));
      m4->getObjectPtr()->csm_host_ = boost::shared_ptr< hoNDArray<float_complext> >( new hoNDArray<float_complext>(*csm_host_));
      m4->getObjectPtr()->reg_host_ = boost::shared_ptr< hoNDArray<float_complext> >( new hoNDArray<float_complext>(*reg_host_));

      GadgetContainerMessage<ISMRMRD::ImageHeader>* m3 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
      m3->cont(m4);

      ISMRMRD::AcquisitionHeader base_head = *m1->getObjectPtr();
      m3->getObjectPtr()->matrix_size[0] = image_dimensions_recon_[0];
      m3->getObjectPtr()->matrix_size[1] = image_dimensions_recon_[1];
      m3->getObjectPtr()->matrix_size[2] = num_frames;
      m3->getObjectPtr()->channels       = num_coils;
      m3->getObjectPtr()->slice          = base_head.idx.slice;
      m3->getObjectPtr()->set            = base_head.idx.set;

      memcpy(m3->getObjectPtr()->position,base_head.position, sizeof(float)*3);
      memcpy(m3->getObjectPtr()->read_dir,base_head.read_dir, sizeof(float)*3);
      memcpy(m3->getObjectPtr()->phase_dir,base_head.phase_dir, sizeof(float)*3);
      memcpy(m3->getObjectPtr()->slice_dir,base_head.slice_dir, sizeof(float)*3);
      memcpy(m3->getObjectPtr()->patient_table_position, base_head.patient_table_position, sizeof(float)*3);

      m3->getObjectPtr()->image_data_type = ISMRMRD::DATA_COMPLEX_FLOAT;
      m3->getObjectPtr()->image_index = ++image_counter_; 
      m3->getObjectPtr()->image_series_index = image_series_;

      //GADGET_DEBUG1("Putting job on queue\n");
      
      if (this->next()->putq(m3) < 0) {
	GADGET_DEBUG1("Failed to put job on queue.\n");
	m3->release();
	return GADGET_FAIL;
      }
    }
    
    if( is_last_profile_in_frame )
      profiles_counter_frame_[set*slices_+slice] = 0;
    else{
      profiles_counter_frame_[set*slices_+slice]++;
    }

    if( new_frame_detected ){

      // This is the first profile of the next frame, enqueue.
      // We have encountered deadlocks if the same profile is enqueued twice in different queues. Hence the copy.
      
      frame_queue_[set*slices_+slice].enqueue_tail(duplicate_profile(m1));
      recon_queue_[set*slices_+slice].enqueue_tail(duplicate_profile(m1)); 

      profiles_counter_frame_[set*slices_+slice]++;
    }

    profiles_counter_global_[set*slices_+slice]++;
    
    m1->release(); // the internal queues hold copies
    return GADGET_OK;
  }
  
  int 
  gpuRadialSenseGadget::calculate_trajectory_for_reconstruction(long profile_offset)
  {   
    //GADGET_DEBUG1("Calculating trajectory for reconstruction\n");

    switch(mode_){
      
    case 0:
    case 1:
      {
	if( rotations_per_reconstruction_ == 0 ){

	  long local_frame = (profile_offset/profiles_per_frame_)%frames_per_rotation_;
	  float angular_offset = M_PI/float(profiles_per_frame_)*float(local_frame)/float(frames_per_rotation_);	  

	  host_traj_recon_ = compute_radial_trajectory_fixed_angle_2d<float>
	    ( samples_per_profile_, profiles_per_frame_, 1, angular_offset )->to_host();	
	}
	else{
	  host_traj_recon_ = compute_radial_trajectory_fixed_angle_2d<float>
	    ( samples_per_profile_, profiles_per_frame_, frames_per_rotation_ )->to_host();
	}
      }
      break;
      
    case 2:
      {
	if( rotations_per_reconstruction_ == 0 ){	  
	  unsigned int first_profile_in_reconstruction = std::max(0L, profile_offset-profiles_per_frame_+1);
	  host_traj_recon_ = compute_radial_trajectory_golden_ratio_2d<float>
	    ( samples_per_profile_, profiles_per_frame_, 1, first_profile_in_reconstruction )->to_host();	
	}
	else{
	  unsigned int first_profile_in_reconstruction = 
	    std::max(0L, profile_offset-profiles_per_frame_*frames_per_rotation_*rotations_per_reconstruction_+1);
	  host_traj_recon_ = compute_radial_trajectory_golden_ratio_2d<float>
	    ( samples_per_profile_, profiles_per_frame_, frames_per_rotation_*rotations_per_reconstruction_, first_profile_in_reconstruction )->to_host();
	}	  
      }
      break;
	
    default:
      GADGET_DEBUG1("Illegal trajectory mode\n");
      return GADGET_FAIL;
      break;
    }
    return GADGET_OK;
  }  

  int
  gpuRadialSenseGadget::calculate_density_compensation_for_reconstruction()
  {
    //GADGET_DEBUG1("Calculating dcw for reconstruction\n");
    
    switch(mode_){
      
    case 0:
    case 1:
      host_weights_recon_ = compute_radial_dcw_fixed_angle_2d<float>
	( samples_per_profile_, profiles_per_frame_, oversampling_factor_, 
	  1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])) )->to_host();
      break;
      
    case 2:
      host_weights_recon_ = compute_radial_dcw_golden_ratio_2d<float>
	( samples_per_profile_, profiles_per_frame_, oversampling_factor_, 
	  1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])) )->to_host();
      break;
      
    default:
      GADGET_DEBUG1("Illegal dcw mode\n");
      return GADGET_FAIL;
      break;
    }
    return GADGET_OK;
  }
  
  boost::shared_ptr< cuNDArray<floatd2> > 
  gpuRadialSenseGadget::calculate_trajectory_for_frame(long profile_offset)
  {
    //GADGET_DEBUG1("Calculating trajectory for buffer frame\n");

    boost::shared_ptr< cuNDArray<floatd2> > result;

    switch(mode_){

    case 0:
    case 1:
      {
	long local_frame = (profile_offset/profiles_per_frame_)%frames_per_rotation_;
	float angular_offset = M_PI/float(profiles_per_frame_)*float(local_frame)/float(frames_per_rotation_);	  

	result = compute_radial_trajectory_fixed_angle_2d<float>
	  ( samples_per_profile_, profiles_per_frame_, 1, angular_offset );  
      }
      break;
	
    case 2:
      { 
	unsigned int first_profile_in_buffer = std::max(0L, profile_offset-profiles_per_frame_+1);

	result  = compute_radial_trajectory_golden_ratio_2d<float>
	  ( samples_per_profile_, profiles_per_frame_, 1, first_profile_in_buffer );
      }
      break;	
	
    default:
      GADGET_DEBUG1("Illegal trajectory mode\n");
      break;
    }
    
    return result;
  }

  boost::shared_ptr< cuNDArray<float> >
  gpuRadialSenseGadget::calculate_density_compensation_for_frame()
  {    
    //GADGET_DEBUG1("Calculating dcw for buffer frame\n");

    switch(mode_){
      
    case 0:
    case 1:
      return compute_radial_dcw_fixed_angle_2d<float>
	( samples_per_profile_, profiles_per_frame_, oversampling_factor_, 1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])) );
      break;
      
    case 2:
      return compute_radial_dcw_golden_ratio_2d<float>
	( samples_per_profile_, profiles_per_frame_, oversampling_factor_, 1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])) );
      break;
      
    default:
      GADGET_DEBUG1("Illegal dcw mode\n");
      return boost::shared_ptr< cuNDArray<float> >();
      break;
    }   
  }


  boost::shared_ptr< cuNDArray<floatd2> > 
  gpuRadialSenseGadget::calculate_trajectory_for_rhs(long profile_offset)
  {
    //GADGET_DEBUG1("Calculating trajectory for rhs\n");

    boost::shared_ptr< cuNDArray<floatd2> > result;

    switch(mode_){

    case 0:
    case 1:
      return compute_radial_trajectory_fixed_angle_2d<float>( samples_per_profile_, profiles_per_frame_*buffer_frames_per_rotation_, 1 );
      break;
	
    case 2:
      { 
	unsigned int first_profile = std::max(0L, profile_offset-profiles_per_frame_*buffer_frames_per_rotation_*buffer_length_in_rotations_+1);
	return compute_radial_trajectory_golden_ratio_2d<float>
	  ( samples_per_profile_, profiles_per_frame_*buffer_frames_per_rotation_*buffer_length_in_rotations_, 1, first_profile );
      }
      break;	
	
    default:
      GADGET_DEBUG1("Illegal trajectory mode\n");
      return boost::shared_ptr< cuNDArray<floatd2> >();
      break;
    }
  }
  
  boost::shared_ptr< cuNDArray<float> >
  gpuRadialSenseGadget::calculate_density_compensation_for_rhs()
  {
    //GADGET_DEBUG1("Calculating dcw for rhs\n");
    
    switch(mode_){
      
    case 0:
    case 1:
      {
	unsigned int num_profiles = profiles_per_frame_*buffer_frames_per_rotation_;
	return compute_radial_dcw_fixed_angle_2d<float>
	  ( samples_per_profile_, num_profiles, oversampling_factor_, 1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])) );
      }
      break;
      
    case 2:
      {
	unsigned int num_profiles = profiles_per_frame_*buffer_frames_per_rotation_*buffer_length_in_rotations_;
	return  compute_radial_dcw_golden_ratio_2d<float>
	  ( samples_per_profile_, num_profiles, oversampling_factor_, 
	    1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])) );
      }
      break;
      
    default:
      GADGET_DEBUG1("Illegal dcw mode\n");
      return boost::shared_ptr< cuNDArray<float> >();
      break;
    }
  }

  boost::shared_ptr< hoNDArray<float_complext> > gpuRadialSenseGadget::extract_samples_from_queue
  ( ACE_Message_Queue<ACE_MT_SYNCH> *queue, std::vector<unsigned int> dims )
  {    
    //GADGET_DEBUG1("Emptying queue...\n");
    unsigned int profiles_buffered = queue->message_count();
    unsigned num_coils = dims.back();
    
    boost::shared_ptr< hoNDArray<float_complext> > host_samples(new hoNDArray<float_complext>(&dims));
    
    for (unsigned int p=0; p<profiles_buffered; p++) {

      ACE_Message_Block* mbq;
      if (queue->dequeue_head(mbq) < 0) {
	GADGET_DEBUG1("Message dequeue failed\n");
	return boost::shared_ptr< hoNDArray<float_complext> >();
      }
      
      GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *acq = AsContainerMessage<ISMRMRD::AcquisitionHeader>(mbq);
      GadgetContainerMessage< hoNDArray< std::complex<float> > > *daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq->cont());
	
      if (!acq || !daq) {
	GADGET_DEBUG1("Unable to interpret data on message queue\n");
	return boost::shared_ptr< hoNDArray<float_complext> >();
      }
	
      for (unsigned int c = 0; c < num_coils; c++) {
	
	float_complext *data_ptr = host_samples->get_data_ptr();
	data_ptr += c*samples_per_profile_*profiles_buffered+p*samples_per_profile_;
	    
	std::complex<float> *r_ptr = daq->getObjectPtr()->get_data_ptr();
	r_ptr += c*daq->getObjectPtr()->get_size(0);
	  
	memcpy(data_ptr,r_ptr,samples_per_profile_*sizeof(float_complext));
      }

      mbq->release();
    }    

    if( queue->message_count() ) {
      GADGET_DEBUG1("Error emptying queue\n");
      return boost::shared_ptr< hoNDArray<float_complext> >();
    }
    
    return host_samples;
  }
  
  GadgetContainerMessage<ISMRMRD::AcquisitionHeader>*
  gpuRadialSenseGadget::duplicate_profile( GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *profile )
  {
    GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *copy = 
      new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();
    
    GadgetContainerMessage< hoNDArray< std::complex<float> > > *cont_copy = 
      new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
    
    *copy->getObjectPtr() = *profile->getObjectPtr();
    *(cont_copy->getObjectPtr()) = *(AsContainerMessage<hoNDArray< std::complex<float> > >(profile->cont())->getObjectPtr());
    
    copy->cont(cont_copy);
    return copy;
  }

  void gpuRadialSenseGadget::reconfigure( unsigned int num_coils )
  {    
    GADGET_DEBUG2("\nReconfiguring:\n#profiles/frame:%d\n#frames/rotation: %d\n#rotations/reconstruction:%d\n", profiles_per_frame_, frames_per_rotation_, rotations_per_reconstruction_);

    calculate_trajectory_for_reconstruction(0);
    calculate_density_compensation_for_reconstruction();
    
    buffer_frames_per_rotation_ = get_int_value(std::string("buffer_frames_per_rotation").c_str());    

    if( buffer_frames_per_rotation_ == 0 ){
      if( mode_ == 2 )
	buffer_frames_per_rotation_ = image_dimensions_recon_os_[0]/profiles_per_frame_;
      else
	buffer_frames_per_rotation_ = frames_per_rotation_;
    }
      
    acc_buffer_->setup( from_std_vector<unsigned int,2>(image_dimensions_recon_), image_dimensions_recon_os_, 
			kernel_width_, num_coils, buffer_length_in_rotations_, buffer_frames_per_rotation_ );
    
    boost::shared_ptr< cuNDArray<float> > device_weights_frame = calculate_density_compensation_for_frame();
    acc_buffer_->set_dcw(device_weights_frame);

    if( buffer_using_solver_ ){
      cuSenseBufferCg<float,2> *acc = (cuSenseBufferCg<float,2>*) acc_buffer_.get();
      acc->set_dcw_for_rhs(calculate_density_compensation_for_rhs());
      acc->preprocess(calculate_trajectory_for_rhs(0).get());
    }
    
    reconfigure_ = false;
  }

  GADGET_FACTORY_DECLARE(gpuRadialSenseGadget)
}
