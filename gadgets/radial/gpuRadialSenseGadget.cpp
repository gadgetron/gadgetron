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

//#include "hoNDArray_elemwise.h"
//#include "hoNDArray_fileio.h"

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
    , host_data_buffer_(0x0)
  {
    GADGET_DEBUG1("Initializing radial base class\n");
    
    // Set some default values in case the config does not contain a specification
    //

    set_parameter(std::string("mode").c_str(), "0");
    set_parameter(std::string("deviceno").c_str(), "0");
    set_parameter(std::string("convolution_kernel_width").c_str(), "5.5");
    set_parameter(std::string("convolution_oversampling_factor").c_str(), "1.25");
    set_parameter(std::string("rotations_per_reconstruction").c_str(), "0");
    set_parameter(std::string("buffer_length_in_rotations").c_str(), "4");
    set_parameter(std::string("real_time_buffer_mode").c_str(), "false");
  }
  
  gpuRadialSenseGadget::~gpuRadialSenseGadget()
  {
    if (host_data_buffer_) delete[] host_data_buffer_;
    if( previous_profile_ ) delete[] previous_profile_;
    if( profiles_counter_frame_ ) delete[] profiles_counter_frame_;
    if( profiles_counter_global_ ) delete[] profiles_counter_global_;
  }

  int gpuRadialSenseGadget::process_config(ACE_Message_Block* mb)
  {
    GADGET_DEBUG1("gpuRadialSenseGadget::process_config\n");
    GPUTimer timer("gpuRadialSenseGadget::process_config");

    mode_ = get_int_value(std::string("mode").c_str());
    device_number_ = get_int_value(std::string("deviceno").c_str());
    profiles_per_frame_ = get_int_value(std::string("profiles_per_frame").c_str());
    frames_per_rotation_ = get_int_value(std::string("frames_per_rotation").c_str());
    rotations_per_reconstruction_ = get_int_value(std::string("rotations_per_reconstruction").c_str());
    buffer_length_in_rotations_ = get_int_value(std::string("buffer_length_in_rotations").c_str());
    cg_iterations_for_reg_ = get_int_value(std::string("cg_iterations_for_regularization").c_str());
    real_time_buffer_mode_ = get_bool_value(std::string("real_time_buffer_mode").c_str());

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

    image_dimensions_.push_back(e_space.matrixSize().x());
    image_dimensions_.push_back(e_space.matrixSize().y());

    if( get_int_value(std::string("reconstruction_matrix_size_x").c_str() ) == 0 )
      image_dimensions_recon_.push_back(e_space.matrixSize().x());      
    else
      image_dimensions_recon_.push_back(get_int_value(std::string("reconstruction_matrix_size_x").c_str()));
    
    if( get_int_value(std::string("reconstruction_matrix_size_y").c_str() ) == 0 )
      image_dimensions_recon_.push_back(e_space.matrixSize().y());      
    else
      image_dimensions_recon_.push_back(get_int_value(std::string("reconstruction_matrix_size_y").c_str()));
    
    GADGET_DEBUG2("matrix_size_x : %d, recon: %d\n", image_dimensions_[0], image_dimensions_recon_[1]);
    GADGET_DEBUG2("matrix_size_y : %d, recon: %d\n", image_dimensions_[1], image_dimensions_recon_[1]);

    // Convolution kernel width and oversampling ratio
    kernel_width_ = get_double_value(std::string("convolution_kernel_width").c_str());
    oversampling_factor_ = get_double_value(std::string("convolution_oversampling_factor").c_str());

    cudaDeviceProp deviceProp;
    if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
      GADGET_DEBUG1( "Error: unable to query device properties.\n" );
      return GADGET_FAIL;
    }
    
    unsigned int warp_size = deviceProp.warpSize;

    // Oversampled matrix size (as a multiple of the warp size, a requirement for the NFFT)
    uintd2 matrix_size = uintd2(image_dimensions_recon_[0],image_dimensions_recon_[1]);
    uintd2 matrix_size_os
      (static_cast<unsigned int>(std::ceil((matrix_size.vec[0]*oversampling_factor_)/warp_size)*warp_size),
       static_cast<unsigned int>(std::ceil((matrix_size.vec[1]*oversampling_factor_)/warp_size)*warp_size));

    oversampling_factor_ = float(matrix_size_os.vec[0])/float(matrix_size.vec[0]); // in case the warp_size constraint above kicked in
    
    // Initialize plan
    plan_ = cuNFFT_plan<float, 2>( matrix_size, matrix_size_os, kernel_width_ );
       
    slices_ = e_limits.slice().present() ? e_limits.slice().get().maximum() + 1 : 1;
    sets_ = e_limits.set().present() ? e_limits.set().get().maximum() + 1 : 1;
    
    previous_profile_ = new long[slices_*sets_];
    profiles_counter_frame_= new long[slices_*sets_];
    profiles_counter_global_= new long[slices_*sets_];

    if( !previous_profile_ || !profiles_counter_frame_ || !profiles_counter_global_ ){
      GADGET_DEBUG1("Failed to allocate host memory");
      return GADGET_FAIL;
    }

    for( unsigned int i=0; i<slices_*sets_; i++ ){
      previous_profile_[i] = -1;
      profiles_counter_frame_[i] = 0;
      profiles_counter_global_[i] = 0;
    }

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

    // Have the imaging plane changed?
    //

    if ( !read_dir_equal(m1->getObjectPtr()->read_dir) || 
	 !phase_dir_equal(m1->getObjectPtr()->phase_dir) ||
	 !slice_dir_equal(m1->getObjectPtr()->slice_dir) ||
	 !position_equal(m1->getObjectPtr()->position)) {

      if (host_data_buffer_){
	delete[] host_data_buffer_;
	host_data_buffer_ = 0x0;
      }

      if( real_time_buffer_mode_ )
	rhs_buffer_->clear();

      memcpy(position_,m1->getObjectPtr()->position,3*sizeof(float));
      memcpy(read_dir_,m1->getObjectPtr()->read_dir,3*sizeof(float));
      memcpy(phase_dir_,m1->getObjectPtr()->phase_dir,3*sizeof(float));
      memcpy(slice_dir_,m1->getObjectPtr()->slice_dir,3*sizeof(float));
    }

    // Allocate host data buffer at the first invocation 
    // - only at this point do we know the actual number of channels after channel reduction

    if (!host_data_buffer_) {
      
      // Assign some default values ("upper bound estimates") of the yet unknown entities. 

      samples_per_profile_ = m1->getObjectPtr()->number_of_samples;
    
      if( profiles_per_frame_ == 0 ){
	profiles_per_frame_ = image_dimensions_[0];
      }

      if( frames_per_rotation_ == 0 ){
	if( mode_ == 2 ) // golden angle
	  frames_per_rotation_ = 1;
	else
	  frames_per_rotation_ = image_dimensions_[0]/profiles_per_frame_;
      }

      // Allocate profile buffer with sufficient storage capabilities
      size_t bsize = (sizeof(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>)+sizeof(GadgetContainerMessage< hoNDArray< std::complex<float> > >))
	*profiles_per_frame_*frames_per_rotation_*(rotations_per_reconstruction_+1)*10;
      
      buffer_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);
      for( unsigned int i=0; i<slices_*sets_; i++ ){
	buffer_[i].high_water_mark(bsize);
	buffer_[i].low_water_mark(bsize);
      }      

      // Allocate memory for profiles buffer (for csm and regularization estimation)
      
      if( prepare_host_buffers(m1->getObjectPtr()->active_channels) == GADGET_FAIL ){
	GADGET_DEBUG1("Failed to create host buffers");
	return GADGET_FAIL;
      }

      if( initialize_regularization_image_solver(m1->getObjectPtr()->active_channels) != GADGET_OK ){
	GADGET_DEBUG1("Failed to initialize regularization solver");
	return GADGET_FAIL;
      }
    }

    unsigned int profile = m1->getObjectPtr()->idx.kspace_encode_step_1;
    unsigned int slice = m1->getObjectPtr()->idx.slice;
    unsigned int set = m1->getObjectPtr()->idx.set;
    bool new_frame_detected = false;

    // Keep track of the incoming profile ids

    if (previous_profile_[set*slices_+slice] >= 0) {

      if ( profile > previous_profile_[set*slices_+slice]) { // this is not the last profile in the frame
	if( mode_ == 0 && get_int_value(std::string("frames_per_rotation").c_str()) == 0 ){
	  unsigned int acceleration_factor = profile - previous_profile_[set*slices_+slice];
	  if( acceleration_factor != frames_per_rotation_ ){
	    frames_per_rotation_ = acceleration_factor;
	    if( prepare_host_buffers(m1->getObjectPtr()->active_channels) == GADGET_FAIL ){
	      GADGET_DEBUG1("Failed to create host buffers");
	      return GADGET_FAIL;
	    }
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
	  if( prepare_host_buffers(m1->getObjectPtr()->active_channels) == GADGET_FAIL ){
	    GADGET_DEBUG1("Failed to create host buffers");
	    return GADGET_FAIL;
	  }
	}
      }
    }
    previous_profile_[set*slices_+slice] = profile;

    if( !real_time_buffer_mode_ ){
    
      // Copy the incoming profile to the host buffer
      //
      
      unsigned int samples_to_copy = m1->getObjectPtr()->number_of_samples;
      unsigned int samples_per_channel =  host_data_buffer_->get_size(0);
      unsigned int profiles_in_buffer = profiles_per_frame_*frames_per_rotation_*buffer_length_in_rotations_;
      
      std::complex<float>* data_ptr = reinterpret_cast< std::complex<float>* >
	(host_data_buffer_[set*slices_+slice].get_data_ptr());
      
      std::complex<float>* profile_ptr = m2->getObjectPtr()->get_data_ptr();
      
      for (unsigned int c = 0; c < m1->getObjectPtr()->active_channels; c++) {
	memcpy(data_ptr+
	       c*samples_per_channel+
	       (profiles_counter_global_[set*slices_+slice]%profiles_in_buffer)*samples_to_copy,
	       profile_ptr+c*m1->getObjectPtr()->number_of_samples, 
	       samples_to_copy*sizeof(std::complex<float>));
      }
    }

    // Are we ready to reconstruct (downstream) ?
    //

    long profiles_per_reconstruction = profiles_per_frame_*frames_per_rotation_*rotations_per_reconstruction_;

    bool is_last_profile_in_frame = (profiles_counter_frame_[set*slices_+slice] == profiles_per_frame_-1);
    is_last_profile_in_frame |= new_frame_detected;
       
    bool is_last_rotation = (rotations_per_reconstruction_ == 0);
    if( profiles_per_reconstruction > 0 )
      is_last_rotation |= ((profiles_counter_global_[set*slices_+slice]%profiles_per_reconstruction) == profiles_per_reconstruction-1);

    // Enqueue profile until we have data for all frames to pass downstream
    //

    if( !new_frame_detected ) // if 'new_frame_detected' we must reconstruct before enqueing the current profile 
      buffer_[set*slices_+slice].enqueue_tail(m1);

    if( is_last_profile_in_frame && is_last_rotation ){
      
      GPUTimer timer("gpuRadialSenseGadget::processing frame(s)");

      //
      // Compute regularization image.
      // -----------------------------
      // This is done either using the cuSenseRHSBuffer (for real time mode)
      // or a host buffer of profiles.
      //

      boost::shared_ptr< cuNDArray<float_complext> > image(new cuNDArray<float_complext>); 

      if( !real_time_buffer_mode_ ){

	// Using the host buffer to compute the regularization image
	//

	if( mode_ == 2 ){
	  
	  // If this is a golden ratio scan, the buffer trajectories (and dcw buffer "permutation") keep changing...
	  //
	  
	  shift_density_compensation_for_buffer(profiles_counter_global_[set*slices_+slice]);
	  
	  boost::shared_ptr< cuNDArray<floatd2> > traj = 
	    calculate_trajectory_for_buffer(profiles_counter_global_[set*slices_+slice]);
	  
	  // Re-preprocess plan
	  try{ plan_.preprocess( traj.get(), cuNFFT_plan<float,2>::NFFT_PREP_ALL ); }
	  catch (runtime_error& err){
	  GADGET_DEBUG_EXCEPTION(err,"NFFT preprocess failed\n");
	  return GADGET_FAIL;
	  }
	}
	
	unsigned int num_coils = m1->getObjectPtr()->active_channels;
	cuNDArray<float_complext> data(&host_data_buffer_[set*slices_+slice]);
	
	// Setup averaged image (an image for each coil)
	std::vector<unsigned int> image_dims;
	image_dims.push_back(image_dimensions_recon_[0]);
	image_dims.push_back(image_dimensions_recon_[1]);
	if( buffer_length_in_rotations_ > 1 && (mode_ == 0 || mode_ == 1) )
	  image_dims.push_back(buffer_length_in_rotations_);
	image_dims.push_back(num_coils);
	
	try{image->create(&image_dims);}
	catch (runtime_error &err){
	  GADGET_DEBUG_EXCEPTION(err,"\nError allocating coil images on device\n");
	  return GADGET_FAIL;
	}
	
	try {
	  plan_.compute( &data, image.get(), device_weights_buffer_.get(), cuNFFT_plan<float,2>::NFFT_BACKWARDS_NC2C );
	}
	catch (runtime_error& err){
	  GADGET_DEBUG_EXCEPTION(err,"NFFT compute failed\n");
	  return GADGET_FAIL;
	}
      	
	// For modes 0 and 1 we have reconstructed coil images for each buffer rotation. Average those.
	if( buffer_length_in_rotations_ > 1 && (mode_ == 0 || mode_ == 1) ){
	  boost::shared_ptr< cuNDArray<float_complext> > tmp = sum( image.get(), 2); 
	  *tmp /= float_complext(buffer_length_in_rotations_);
	  image = tmp;	
	}
      }
      else{
	
	// Real-time mode using the cuSenseRHSBuffer
	//

	
      }

      // Estimate CSM
      boost::shared_ptr< cuNDArray<float_complext> > csm  = estimate_b1_map<float,2>( image.get() );
      
      // Compute regularization image
      boost::shared_ptr< cuNDArray<float_complext> > reg_image = 
	compute_regularization_image( image.get(), &data, csm, set, slice );
      
      boost::shared_ptr< hoNDArray<float_complext> > csm_host = csm->to_host();
      boost::shared_ptr< hoNDArray<float_complext> > reg_host = reg_image->to_host();
      
      csm.reset();
      reg_image.reset();

      /*
      static int counter = 0;
      char filename[256];
      sprintf((char*)filename, "reg_%d.real", counter);
      write_nd_array<float>( abs(reg_host.get()).get(), filename );
      counter++; */
      
      // Prepare data array of the profiles for the downstream reconstruction
      //
      
      unsigned int profiles_buffered = buffer_[set*slices_+slice].message_count();
      unsigned int num_frames = std::max(1L,frames_per_rotation_*rotations_per_reconstruction_);

      std::vector<unsigned int> ddimensions;
      ddimensions.push_back(samples_per_profile_*profiles_per_frame_*num_frames);
      ddimensions.push_back(num_coils);
      
      boost::shared_ptr< hoNDArray<float_complext> > data_host(new hoNDArray<float_complext>());
      
      try{data_host->create(&ddimensions);}
      catch (runtime_error& err){
	GADGET_DEBUG_EXCEPTION(err,"Unable to allocate host data array\n");
	return GADGET_FAIL;
      }
      
      ddimensions.clear();
      ddimensions.push_back(samples_per_profile_*profiles_per_frame_);
      ddimensions.push_back(num_frames);
      
      for (unsigned int p = 0; p < profiles_buffered; p++) {
	ACE_Message_Block* mbq;
	if (buffer_[set*slices_+slice].dequeue_head(mbq) < 0) {
	  GADGET_DEBUG1("Message dequeue failed\n");
	  return GADGET_FAIL;
	}
	
	GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* acq = AsContainerMessage<ISMRMRD::AcquisitionHeader>(mbq);
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq->cont());
	  
	if (!acq || !daq) {
	  GADGET_DEBUG1("Unable to interpret data on message Q\n");
	  return GADGET_FAIL;
	}
	  
	for (unsigned int c = 0; c < num_coils; c++) {

	  float_complext* data_ptr = data_host->get_data_ptr();
	  data_ptr += c*samples_per_profile_*profiles_buffered+p*samples_per_profile_;
	    
	  std::complex<float>* r_ptr = daq->getObjectPtr()->get_data_ptr();
	  r_ptr += c*daq->getObjectPtr()->get_size(0);
	    
	  memcpy(data_ptr,r_ptr,samples_per_profile_*sizeof(float_complext));
	}
	  
	mbq->release();
      }
	
      if (buffer_[set*slices_+slice].message_count()) {
	GADGET_DEBUG1("Error occured, all messages should have been cleared off the buffer by now.\n");
	return GADGET_FAIL;
      }

      // The trajectory needs to be updated on the fly:
      // - for golden ratio acquisitions
      // - when we are reconstructing frame-by-frame
      
      if( mode_ == 2 || rotations_per_reconstruction_ == 0 ){
	calculate_trajectory_for_reconstruction( profiles_counter_global_[set*slices_+slice] - ((new_frame_detected) ? 1 : 0) );
      }
      
      GadgetContainerMessage< SenseJob >* m4 = new GadgetContainerMessage< SenseJob >();
	
      m4->getObjectPtr()->dat_host_ = data_host;
      m4->getObjectPtr()->csm_host_ = csm_host;
      m4->getObjectPtr()->reg_host_ = reg_host;
      m4->getObjectPtr()->tra_host_ = boost::shared_ptr< hoNDArray<floatd2> >(new hoNDArray<floatd2>(*host_traj_frame_));
      m4->getObjectPtr()->dcw_host_ = boost::shared_ptr< hoNDArray<float> >(new hoNDArray<float>(*host_weights_frame_));

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

      GADGET_DEBUG1("Putting job on queue\n");
      
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
      buffer_[set*slices_+slice].enqueue_tail(m1); // This is the first profile of the next frame
      profiles_counter_frame_[set*slices_+slice]++;
    }

    profiles_counter_global_[set*slices_+slice]++;
    
    return GADGET_OK;
  }

  int gpuRadialSenseGadget::prepare_host_buffers(unsigned int num_channels)
  {    
    GADGET_DEBUG1("Allocating host buffers\n");
    GADGET_DEBUG2("\nsamples_per_profile: %d\nprofiles_per_frame: %d\nframes_per_rotation: %d\nbuffer_length_in_rotations: %d\n", 
		  samples_per_profile_, profiles_per_frame_, frames_per_rotation_, buffer_length_in_rotations_);

    try{
      
      unsigned int elements_per_channel = samples_per_profile_*profiles_per_frame_*frames_per_rotation_*buffer_length_in_rotations_;

      std::vector<unsigned int> data_dimensions;
      data_dimensions.push_back(elements_per_channel);
      data_dimensions.push_back(num_channels);

      // Keep a copy in order to copy over the existing profiles to the new buffer
      hoNDArray<float_complext> *host_data_buffer_old = host_data_buffer_;
      host_data_buffer_ = new hoNDArray<float_complext>[slices_*sets_];
      
      if (!host_data_buffer_) {
	GADGET_DEBUG1("Unable to allocate array for host data buffer\n");
	if( host_data_buffer_old ) delete[] host_data_buffer_old;
	return GADGET_FAIL;
      }
      
      for (unsigned int i = 0; i < slices_*sets_; i++) {
	host_data_buffer_[i].create(&data_dimensions);
	host_data_buffer_[i].fill(0.0f);
      }       
      
      // Copy over profiles from the old buffer
      if( host_data_buffer_old ){
	for (unsigned int i = 0; i < slices_*sets_; i++) {
	  for( unsigned int c=0; c<num_channels; c++ ){
	    memcpy( host_data_buffer_[i].get_data_ptr() + c*elements_per_channel,
		    host_data_buffer_old[i].get_data_ptr() + c*host_data_buffer_old[i].get_size(0),
		    std::min(elements_per_channel, host_data_buffer_old[i].get_size(0))*sizeof(float_complext) );
	  }
	}
	delete[] host_data_buffer_old;
      }
      
      // Initialise density compensation weights
      device_weights_buffer_unpermuted_.reset();
      calculate_density_compensation_for_buffer();
      calculate_density_compensation_for_reconstruction();
      
      // Initialize trajectories
      calculate_trajectory_for_reconstruction(0);
      boost::shared_ptr< cuNDArray<floatd2> > traj = calculate_trajectory_for_buffer(0);
      
      // Preprocess plan (now that we have the trajectory)
      plan_.preprocess( traj.get(), cuNFFT_plan<float,2>::NFFT_PREP_ALL );

      for( unsigned int i=0; i<slices_*sets_; i++ ){
	profiles_counter_global_[i] = profiles_counter_frame_[i];
      }
    }
    
    catch (runtime_error& err){
      GADGET_DEBUG_EXCEPTION(err,"gpuRadialSenseGadget::prepare_host_buffers() failed\n");
      return GADGET_FAIL;
    }

    return GADGET_OK;
  }

  boost::shared_ptr< cuNDArray<floatd2> > 
  gpuRadialSenseGadget::calculate_trajectory_for_buffer(long profile_offset)
  {
    GADGET_DEBUG1("Calculating trajectory for buffer\n");

    boost::shared_ptr< cuNDArray<floatd2> > result;

    switch(mode_){

      case 0:
      case 1:
	{
	  result = compute_radial_trajectory_fixed_angle_2d<float>( samples_per_profile_, profiles_per_frame_, frames_per_rotation_ );
	  
	  // Expand to the number of buffer rotations
	  std::vector<unsigned int> dims;
	  dims.push_back(samples_per_profile_*profiles_per_frame_*frames_per_rotation_);	  
	  result->reshape(&dims);
	  result = expand(result.get(), buffer_length_in_rotations_);
	}
	break;
	
      case 2:
	{ 
	  unsigned int profiles_in_buffer = profiles_per_frame_*frames_per_rotation_*buffer_length_in_rotations_;
	  unsigned int first_profile_in_buffer = std::max(0L, profile_offset-profiles_in_buffer+1);

	  boost::shared_ptr< cuNDArray<floatd2> > tmp = compute_radial_trajectory_golden_ratio_2d<float>
	    ( samples_per_profile_, profiles_in_buffer, 1, first_profile_in_buffer );
	  
	  if( (first_profile_in_buffer%profiles_in_buffer) == 0 )
	    result = tmp;
	  else
	    result = wrap_array( tmp, profile_offset );
	}
	break;	
	
      default:
	GADGET_DEBUG1("Illegal trajectory mode\n");
	break;
      }
        
    return result;
  }

  int 
  gpuRadialSenseGadget::calculate_trajectory_for_reconstruction(long profile_offset)
  {   
    GADGET_DEBUG1("Calculating trajectory for reconstruction\n");

    switch(mode_){
      
    case 0:
    case 1:
      {
	if( rotations_per_reconstruction_ == 0 ){

	  long local_frame = (profile_offset/profiles_per_frame_)%frames_per_rotation_;
	  float angular_offset = M_PI/float(profiles_per_frame_)*float(local_frame)/float(frames_per_rotation_);	  

	  host_traj_frame_ = compute_radial_trajectory_fixed_angle_2d<float>
	    ( samples_per_profile_, profiles_per_frame_, 1, angular_offset )->to_host();	
	}
	else{
	  host_traj_frame_ = compute_radial_trajectory_fixed_angle_2d<float>
	    ( samples_per_profile_, profiles_per_frame_, frames_per_rotation_ )->to_host();
	}
      }
      break;
      
    case 2:
      {
	if( rotations_per_reconstruction_ == 0 ){	  
	  unsigned int first_profile_in_reconstruction = std::max(0L, profile_offset-profiles_per_frame_+1);
	  host_traj_frame_ = compute_radial_trajectory_golden_ratio_2d<float>
	    ( samples_per_profile_, profiles_per_frame_, 1, first_profile_in_reconstruction )->to_host();	
	}
	else{
	  unsigned int first_profile_in_reconstruction = 
	    std::max(0L, profile_offset-profiles_per_frame_*frames_per_rotation_*rotations_per_reconstruction_+1);
	  host_traj_frame_ = compute_radial_trajectory_golden_ratio_2d<float>
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
  gpuRadialSenseGadget::calculate_density_compensation_for_buffer()
  {    
    GADGET_DEBUG1("Calculating dcw for buffer\n");

    switch(mode_){
      
    case 0:
    case 1:
      {
	unsigned int profiles_in_buffer = profiles_per_frame_*frames_per_rotation_;	

	device_weights_buffer_ = compute_radial_dcw_fixed_angle_2d<float>
	  ( samples_per_profile_, profiles_in_buffer, oversampling_factor_, 1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])) );
      }
      break;
      
    case 2:
      {
	unsigned int profiles_in_buffer = profiles_per_frame_*frames_per_rotation_*buffer_length_in_rotations_;
	
	device_weights_buffer_unpermuted_ = compute_radial_dcw_golden_ratio_2d<float>
	  ( samples_per_profile_, profiles_in_buffer, oversampling_factor_, 1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])) );
	
	device_weights_buffer_ = boost::shared_ptr< cuNDArray<float> >(new cuNDArray<float>(device_weights_buffer_unpermuted_.get()));
      }
      break;
      
    default:
      GADGET_DEBUG1("Illegal dcw mode\n");
      return GADGET_FAIL;
    }   
    return GADGET_OK;
  }

  int 
  gpuRadialSenseGadget::shift_density_compensation_for_buffer(long profile_offset)
  {    
    GADGET_DEBUG1("shifting dcw for buffer\n");

    switch(mode_){
      
    case 2:
      {
	if( device_weights_buffer_unpermuted_.get() == 0x0 ){
	  GADGET_DEBUG1("Cannot shift dcw, as it is not yet computed");
	  return GADGET_FAIL;
	}
	
	unsigned int profiles_in_buffer = profiles_per_frame_*frames_per_rotation_*buffer_length_in_rotations_;
	unsigned int first_profile_in_buffer = std::max(0L, profile_offset-profiles_in_buffer+1);
	if( (first_profile_in_buffer%profiles_in_buffer) != 0 )	  
	  device_weights_buffer_ = wrap_array( device_weights_buffer_unpermuted_, profile_offset );
      }
      break;
      
    default:
      GADGET_DEBUG1("Illegal dcw mode\n");
      return GADGET_FAIL;
      break;
    }

    return GADGET_OK;
  }

  int
  gpuRadialSenseGadget::calculate_density_compensation_for_reconstruction()
  {
    GADGET_DEBUG1("Calculating dcw for reconstruction\n");
    
    switch(mode_){
      
    case 0:
    case 1:
      {
	host_weights_frame_ = compute_radial_dcw_fixed_angle_2d<float>
	  ( samples_per_profile_, profiles_per_frame_, oversampling_factor_, 
	    1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])) )->to_host();
      }
      break;
      
    case 2:
      {
	host_weights_frame_ = compute_radial_dcw_golden_ratio_2d<float>
	  ( samples_per_profile_, profiles_per_frame_, oversampling_factor_, 
	    1.0f/(float(samples_per_profile_)/float(image_dimensions_recon_[0])) )->to_host();
      }
      break;
      
    default:
      GADGET_DEBUG1("Illegal dcw mode\n");
      return GADGET_FAIL;
      break;
    }
    return GADGET_OK;
  }

  template<class ARRAY> boost::shared_ptr<ARRAY> 
  gpuRadialSenseGadget::wrap_array( boost::shared_ptr<ARRAY> tmp, long profile_offset )
  {
    GADGET_DEBUG1("shifting host buffer contents\n");

    // The temporary array ('tmp') should be "permuted" to fit the order of the profiles in the buffer:
    // - the first sample in 'tmp' corresponds to the position of buffer profile 'profile_offset' modulus the buffer size
    // - and 'tmp' then "wraps around" the buffer
    
    boost::shared_ptr<ARRAY> result(new ARRAY(tmp->get_dimensions()));

    unsigned int profiles_in_buffer = profiles_per_frame_*frames_per_rotation_*buffer_length_in_rotations_;
    long profile_idx_in_buffer = profile_offset%profiles_in_buffer;
    
    std::vector<unsigned int> head_out_dimensions, tail_out_dimensions;
    head_out_dimensions.push_back((profile_idx_in_buffer+1)*samples_per_profile_);
    tail_out_dimensions.push_back((profiles_in_buffer-profile_idx_in_buffer-1)*samples_per_profile_);
    
    ARRAY head_out( &head_out_dimensions, result->get_data_ptr() );
    ARRAY tail_out( &tail_out_dimensions, result->get_data_ptr()+(profile_idx_in_buffer+1)*samples_per_profile_ );
    
    ARRAY head_in( &tail_out_dimensions, tmp->get_data_ptr() );
    ARRAY tail_in( &head_out_dimensions, tmp->get_data_ptr()+(profiles_in_buffer-profile_idx_in_buffer-1)*samples_per_profile_ );
    
    head_out = tail_in;
    tail_out = head_in;		
    
    return result;
  }
  

  int
  gpuRadialSenseGadget::initialize_regularization_image_solver(unsigned int num_channels)
  {
    // Setup solver
    
    E_ = boost::shared_ptr< cuNonCartesianSenseOperator<float,2> >( new cuNonCartesianSenseOperator<float,2>() );        
    
    if( real_time_buffer_mode_ ){
      rhs_buffer_ = boost::shared_ptr< cuSenseRHSBuffer<float,2> >( new cuSenseRHSBuffer<float,2>() );
      rhs_buffer_->set_num_coils( num_channels );
      rhs_buffer_->set_sense_operator( E_ );
    }

    if( (cg_iterations_for_reg_ > 0) || real_time_buffer_mode_ ){ 
      
      cudaDeviceProp deviceProp;
      if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
	GADGET_DEBUG1( "Error: unable to query device properties.\n" );
	return GADGET_FAIL;
      }
      
      unsigned int warp_size = deviceProp.warpSize;
      
      uintd2 matrix_size = uintd2(image_dimensions_recon_[0],image_dimensions_recon_[1]);
      uintd2 matrix_size_os
	(static_cast<unsigned int>(std::ceil((matrix_size.vec[0]*oversampling_factor_)/warp_size)*warp_size),
	 static_cast<unsigned int>(std::ceil((matrix_size.vec[1]*oversampling_factor_)/warp_size)*warp_size));
      
      try{ E_->setup( matrix_size, matrix_size_os, kernel_width_); }
      catch (runtime_error& err){
	GADGET_DEBUG_EXCEPTION(err, "\nError: unable to setup encoding operator.\n" );
	return GADGET_FAIL;
      }
      
      E_->set_dcw(device_weights_buffer_);      
      
      D_ = boost::shared_ptr< cuCgPreconditioner<float_complext> >( new cuCgPreconditioner<float_complext>() );
      
      cg_.set_encoding_operator( E_ );        // encoding matrix
      cg_.set_preconditioner( D_ );           // preconditioning matrix
      cg_.set_max_iterations( cg_iterations_for_reg_ );
      cg_.set_tc_tolerance( 1e-6 );
      cg_.set_output_mode( cuCgSolver<float_complext>::OUTPUT_SILENT);
    }
    return GADGET_OK;
  }

  boost::shared_ptr< cuNDArray<float_complext> > 
  gpuRadialSenseGadget::compute_regularization_image( cuNDArray<float_complext> *image, cuNDArray<float_complext> *data, 
						      boost::shared_ptr< cuNDArray<float_complext> > csm, 
						      unsigned int set, unsigned int slice )
  {
    boost::shared_ptr< cuNDArray<float_complext> > reg_image;
    
    try{
      
      // Setup averaged image

      std::vector<unsigned int> image_dims;
      image_dims.push_back(image_dimensions_recon_[0]); image_dims.push_back(image_dimensions_recon_[1]);
      if( buffer_length_in_rotations_ > 1 && (mode_ == 0 || mode_ == 1) )
	image_dims.push_back(buffer_length_in_rotations_);
      
      boost::shared_ptr< std::vector<unsigned int> > reg_dims = image->get_dimensions(); reg_dims->pop_back();    
      reg_image = boost::shared_ptr< cuNDArray<float_complext> >( new cuNDArray<float_complext>(reg_dims.get()));
      
      E_->set_csm(csm);

      // There are two options of computing the regularization image:
      //
      
      if( cg_iterations_for_reg_ == 0 ) {
	
	// 1. Basic coil combination based on the csm (the default mode)
	// 
	
	E_->mult_csm_conj_sum( image, reg_image.get() );
      }
      else{
	
	// 2. Compute regularization based on a conjugate gradient solver
	//    - preferred mode if the scaling of the regularization image should match the reconstructed images
	// 
			
	E_->set_domain_dimensions(&image_dims);
	E_->set_codomain_dimensions(data->get_dimensions().get());
	
	if( !E_->is_preprocessed() || mode_ == 2 ){ // mode 2 is golden ratio

	  boost::shared_ptr< cuNDArray<floatd2> > traj = 
	    calculate_trajectory_for_buffer(profiles_counter_global_[set*slices_+slice]);
	  E_->preprocess(traj.get());
	}
	
	// Define preconditioning weights
	boost::shared_ptr< cuNDArray<float> > _precon_weights = sum(abs_square(csm.get()).get(), 2);
	reciprocal_sqrt_inplace(_precon_weights.get());	
	boost::shared_ptr< cuNDArray<float_complext> > precon_weights = real_to_complex<float_complext>( _precon_weights.get() );
	_precon_weights.reset();
	D_->set_weights( precon_weights ); 
	
	// Compute regularization image
	boost::shared_ptr< cuNDArray<float_complext> > cgresult = cg_.solve(data);
	
	// For modes 0 and 1 we have reconstructed coil images for each buffer rotation. Average those.
	if( buffer_length_in_rotations_ > 1 && (mode_ == 0 || mode_ == 1) ){
	  boost::shared_ptr< cuNDArray<float_complext> > tmp = sum( cgresult.get(), 2); 
	  *tmp /= float_complext(buffer_length_in_rotations_);
	  *reg_image = *tmp;	
	}
	else
	  *reg_image = *cgresult;
	
	/*
	static int counter = 0;
	char filename[256];
	sprintf((char*)filename, "cgreg_%d.real", counter);
	write_nd_array<float>( abs(cgresult.get())->to_host().get(), filename );
	counter++; */
      }
    }
    catch (runtime_error& err){
      GADGET_DEBUG_EXCEPTION(err,"Error computing regularization image");
    }
    return reg_image;
  }
    
  GADGET_FACTORY_DECLARE(gpuRadialSenseGadget)
}
