#include "GPUCGGadget.h"
#include "ConfigParser.h"




GPUCGGadget::GPUCGGadget()
  : profiles_per_frame_(64)
  , shared_profiles_(0)
  , channels_(0)
  , samples_per_profile_(0)
  , device_number_(0)
  , number_of_iterations_(10)
  , oversampling_(1.25)
  , kernel_width_(5.5)
  , kappa_(0.1)
  , gc_factor_(1.0)
  , is_configured_(false)
  , current_profile_offset_(0)
  , current_frame_number_(0)
  , allocated_samples_(0)
  , data_host_ptr_(0)
  , data_dev_ptr_(0)
  , trajectory_dev_ptr_(0)
  , dcw_dev_ptr_(0)
  , csm_buffer_dev_ptr_(0)
  , csm_acc_coil_image_os_dev_ptr_(0)
  , kspace_acc_coil_images_os_size_(0)
  , image_dev_ptr_(0)
  , csm_buffer_length_(8)

{
  matrix_size_    = make_uint2(0,0);
  matrix_size_os_ = make_uint2(0,0);
}


GPUCGGadget::~GPUCGGadget()
{  
  if (data_host_ptr_)      delete [] data_host_ptr_;
  if (data_dev_ptr_)       cudaFree(data_dev_ptr_);
  if (trajectory_dev_ptr_) cudaFree(trajectory_dev_ptr_);
  if (dcw_dev_ptr_)        cudaFree(dcw_dev_ptr_);
}

int GPUCGGadget::process_config(ACE_Message_Block* mb)
{
  GADGET_DEBUG1("GPUCGGadget::process_config\n");

  ConfigParser cp;
  cp.parse(mb->rd_ptr());

  if (!is_configured_) {
    //Initialize Cuda
    cudaDeviceProp deviceProp;
    unsigned int device = device_number_;
    cudaGetDeviceProperties(&deviceProp, device);			
    if (deviceProp.major < 1) {					
      GADGET_DEBUG1("GPU Device does not support CUDA\n");
      return GADGET_FAIL;					
    }
    GADGET_DEBUG2("Using GPU device %d, %s\n", device, deviceProp.name);
    cudaSetDevice(device);						
    cublasInit();
    //End of Cuda Initilization

    samples_per_profile_ = cp.getIntVal("encoding","readout_length");
    channels_ = cp.getIntVal("encoding","channels");

    //Only set matrix size if not set from the outside
    if (matrix_size_.x == 0 && matrix_size_.y == 0) {
      matrix_size_ = make_uint2(cp.getIntVal("encoding","matrix_x"), 
				cp.getIntVal("encoding","matrix_y"));
    }

    GADGET_DEBUG2("Matrix size  : [%d,%d] \n", matrix_size_.x, matrix_size_.y);

    matrix_size_os_ = 
      make_uint2(static_cast<unsigned int>(ceil((matrix_size_.x*oversampling_)/32.0f)*32),
		 static_cast<unsigned int>(ceil((matrix_size_.y*oversampling_)/32.0f)*32));

    GADGET_DEBUG2("Matrix size OS: [%d,%d] \n", matrix_size_os_.x, matrix_size_os_.y);
    
    domain_size_grid_ = make_uint2( 1, 1 );
    domain_size_samples_ = 1;
  
    if (channels_%2 && channels_ > 8) {
      GADGET_DEBUG1("Odd number of coils and more than 8 coils detected. Bailing out\n");
      return GADGET_FAIL;
    }

    if (channels_ <= 8) {
      domain_size_coils_ = channels_;
    } else {
      int groups = 1;
      while (((channels_ / groups) > 8) || (channels_%groups)) {
	groups++;
      }
      domain_size_coils_ = channels_/groups;
      GADGET_DEBUG2("Gridding convolution will be split into %d runs (groups)\n", groups);
      GADGET_DEBUG2("      Domain size: %d\n", domain_size_coils_);
    }

    
    fixed_dims_ = make_uint2( 0, 0 );
    
    if (calculate_trajectory() == GADGET_FAIL) {
      GADGET_DEBUG1("Unable to calculate trajectory\n");
      return GADGET_FAIL;
    }

    plan_generic_ = 
      preprocess_generic_NFFT< uint2, float2 >( matrix_size_, 
						matrix_size_os_, 
						fixed_dims_, 
						domain_size_grid_, 
						domain_size_samples_, 
						domain_size_coils_, 
						kernel_width_, 
						profiles_per_frame_*samples_per_profile_, 
						trajectory_dev_ptr_);

    iterative_sense_initialize( plan_generic_, channels_ );

    cudaError_t err = cudaGetLastError();
    if( err != cudaSuccess ){
      GADGET_DEBUG2("Unable to initialize iterative SENSE: %s\n",
		    cudaGetErrorString(err));
    
      
      return GADGET_FAIL;
    }
    
    if (allocate_csm_buffer() == GADGET_FAIL) {
      GADGET_DEBUG1("Unable to allocate CSM buffer\n");
      return GADGET_FAIL;
    }


    if (image_dev_ptr_) cudaFree(image_dev_ptr_);
    cudaMalloc( (void**) &image_dev_ptr_, prod(matrix_size_)*sizeof(cuFloatComplex) );

    err = cudaGetLastError();
    if( err != cudaSuccess ){
      GADGET_DEBUG2("Failed to allocate memory for image: %s\n",
		    cudaGetErrorString(err));
    
      
      return GADGET_FAIL;
    }

    is_configured_ = true;
  }


  return 0;
}


int  GPUCGGadget::process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
			  GadgetContainerMessage< NDArray< std::complex<float> > >* m2)
{
  if (!is_configured_) {
    GADGET_DEBUG1("Data received before configuration complete\n");
    return GADGET_FAIL;
  }

  buffer_.enqueue_tail(m1);
  
  if ((int)buffer_.message_count() >= profiles_per_frame_) {
    //GADGET_DEBUG2("We now have enough data for reconstructing a frame (%d)\n",
    //		  buffer_.message_count());


    if (calculate_trajectory() == GADGET_FAIL) {
      GADGET_DEBUG1("Failed to calculate the trajectory on the GPU\n");
      return GADGET_FAIL;
    }
    
    if (calculate_density_compensation() == GADGET_FAIL) {
      GADGET_DEBUG1("Failed to calculate the density compensation on the GPU\n");
      return GADGET_FAIL;
    }

    if (upload_samples() == GADGET_FAIL) {
      GADGET_DEBUG1("Failed to upload samples to the GPU\n");
      return GADGET_FAIL;
    }

    

    int samples_in_frame = profiles_per_frame_*samples_per_profile_;

    if (!preprocess_generic_NFFT( plan_generic_, 
				  samples_in_frame, 
				  trajectory_dev_ptr_)) {
      GADGET_DEBUG1("Error during preprocess_generic_NFFT\n");
      return GADGET_FAIL;
      
    }


    float shutter_radius = 0.95*0.5;
    if (!noise_decorrelate_generic(samples_in_frame, 
				   channels_, 
				   shutter_radius, 
				   data_dev_ptr_, 
				   trajectory_dev_ptr_ )) {
      
      GADGET_DEBUG1("Error during noise decorrelation\n");
      return GADGET_FAIL;
    } 
    
    
    /*
      if (m_csm_needs_reset) {
      AllocateCSMBuffer();
      }
    */
    
    float sigma_csm = 16.0f;
    unsigned long ptr_offset = 
      (current_frame_number_%csm_buffer_length_)*prod(matrix_size_os_)*channels_;
    
    if (!update_csm_and_regularization( plan_generic_, 
					data_dev_ptr_,
					trajectory_dev_ptr_, 
					dcw_dev_ptr_,
					sigma_csm, 
					&csm_buffer_dev_ptr_[ptr_offset],
					csm_acc_coil_image_os_dev_ptr_,
					kspace_acc_coil_images_os_size_, 
					true,
					false )) { //csm_needs_reset
      
      GADGET_DEBUG1("update_csm_and_regularization failed\n");
      return GADGET_FAIL;    
    }
    
    set_dcw(plan_generic_,dcw_dev_ptr_);
    
    if (!iterative_sense_compute( plan_generic_, 
				  number_of_iterations_, 
				  kappa_, 
				  data_dev_ptr_, 
				  image_dev_ptr_ )) {
      GADGET_DEBUG1("iterative_sense_compute failed\n");
      return GADGET_FAIL;
    }
  

    //Now pass the reconstructed image on
    GadgetContainerMessage<GadgetMessageImage>* cm1 = 
      new GadgetContainerMessage<GadgetMessageImage>();
    
    GadgetContainerMessage< NDArray< std::complex<float> > >* cm2 = 
      new GadgetContainerMessage<NDArray< std::complex<float> > >();
    
    cm1->cont(cm2);
    
    std::vector<int> img_dims(2);
    img_dims[0] = matrix_size_.x;
    img_dims[1] = matrix_size_.y;
    
    if (!cm2->getObjectPtr()->create(img_dims)) {
      GADGET_DEBUG1("Unable to allocate new image array");
      cm1->release();
      return GADGET_FAIL;
    }
    
    size_t data_length = prod(matrix_size_);
    
    cudaMemcpy(cm2->getObjectPtr()->get_data_ptr(),
	       image_dev_ptr_,
	       data_length*sizeof(cuFloatComplex),
	       cudaMemcpyDeviceToHost);

    cm1->getObjectPtr()->matrix_size[0] = img_dims[0];
    cm1->getObjectPtr()->matrix_size[1] = img_dims[1];
    cm1->getObjectPtr()->matrix_size[2] = 1;
    cm1->getObjectPtr()->channels       = 1;
    cm1->getObjectPtr()->data_idx_min       = m1->getObjectPtr()->min_idx;
    cm1->getObjectPtr()->data_idx_max       = m1->getObjectPtr()->max_idx;
    cm1->getObjectPtr()->data_idx_current   = m1->getObjectPtr()->idx;	
    
    memcpy(cm1->getObjectPtr()->position,m1->getObjectPtr()->position,
	   sizeof(float)*3);

    memcpy(cm1->getObjectPtr()->quarternion,m1->getObjectPtr()->quarternion,
	   sizeof(float)*4);
 
    if (this->next()->putq(cm1) < 0) {
      GADGET_DEBUG1("Failed to result image on to Q\n");
      cm1->release();
      return GADGET_FAIL;
    }
    
    //GADGET_DEBUG2("Frame %d reconstructed an passed down the chain\n",
    //current_frame_number_);

    //Dequeue the message we don't need anymore
    ACE_Message_Block* mb_tmp;
    for (int i = 0; i < (profiles_per_frame_-shared_profiles_); i++) {
      buffer_.dequeue_head(mb_tmp);
      mb_tmp->release();
      current_profile_offset_++; 
    }
    
    current_frame_number_++;
  }

  return GADGET_OK;
}

int GPUCGGadget::calculate_trajectory()
{
  if (trajectory_dev_ptr_) {
    cudaFree(trajectory_dev_ptr_);
    trajectory_dev_ptr_ = 0;
  }
  
  if (!trajectory_dev_ptr_) {
    trajectory_dev_ptr_ = compute_trajectory_radial_2d<1>(matrix_size_.x, 
							  matrix_size_os_.x, 
							  samples_per_profile_, 
							  profiles_per_frame_, 
							  current_profile_offset_, 
							  1,
							  gc_factor_);
  }

  if (!trajectory_dev_ptr_) {
    GADGET_DEBUG1("Failed to allocate trajectory");
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

int GPUCGGadget::calculate_density_compensation()
{

  //TODO: add check to see if we really need to recalculate this.
  if (dcw_dev_ptr_) {
    cudaFree(dcw_dev_ptr_);
    dcw_dev_ptr_ = 0;
  }

  if (!dcw_dev_ptr_) {
    dcw_dev_ptr_ = compute_dcw_radial_2d<1>( matrix_size_.x, 
					     matrix_size_os_.x, 
					     samples_per_profile_, 
					     profiles_per_frame_, 
					     current_profile_offset_, 
					     1,
					     gc_factor_);
  }

  if (!dcw_dev_ptr_) {
    GADGET_DEBUG1("Failed to calculate density compensation weights\n");
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

int GPUCGGadget::upload_samples()
{
  int samples_needed = 
    samples_per_profile_*
    profiles_per_frame_;

  if (samples_needed != allocated_samples_) {
    if (data_dev_ptr_) cudaFree(data_dev_ptr_);

    cudaMalloc( (void**) &data_dev_ptr_, 
		channels_*samples_needed*sizeof(cuFloatComplex) );

    cudaError_t err = cudaGetLastError();
    if( err != cudaSuccess ){
      GADGET_DEBUG2("Unable to allocate GPU memory for samples: %s",
		   cudaGetErrorString(err));

      return GADGET_FAIL;
    }

    try {
      data_host_ptr_ = new float[channels_*samples_needed*2];
    } catch (...) {
      GADGET_DEBUG1("Failed to allocate host memory for samples\n");
      return GADGET_FAIL;
    }

    allocated_samples_ = samples_needed;

  }

  ACE_Message_Queue_Reverse_Iterator<ACE_MT_SYNCH> it(buffer_);
  int profiles_copied = 0;
  GadgetContainerMessage<GadgetMessageAcquisition>* m1;
  GadgetContainerMessage< NDArray< std::complex<float> > >* m2;
  ACE_Message_Block* mb;

  while (profiles_copied < profiles_per_frame_) {
    it.next(mb);
    
    m1 = dynamic_cast< GadgetContainerMessage< GadgetMessageAcquisition >* >(mb);
    if (!m1) {
      GADGET_DEBUG1("Failed to dynamic cast message\n");
      return -1;
    }

      
    m2 = dynamic_cast< GadgetContainerMessage< NDArray< std::complex<float> > >* > (m1->cont());
  
    if (!m2) {
      GADGET_DEBUG1("Failed to dynamic cast message\n");
      return -1;
    }
    
    
    std::complex<float>* d = m2->getObjectPtr()->get_data_ptr();
    int current_profile = profiles_per_frame_-profiles_copied-1;
    
    for (int i = 0; i < channels_; i++) {
      
      memcpy(data_host_ptr_ + 
	     (i*allocated_samples_ + current_profile*samples_per_profile_) * 2,
	     d + i*samples_per_profile_, sizeof(float)*samples_per_profile_*2);
      
    }
    
    
    it.advance();   
    profiles_copied++;
  }

  cudaMemcpy( data_dev_ptr_,
	      data_host_ptr_,
	      samples_needed*channels_*sizeof(cuFloatComplex), 
	      cudaMemcpyHostToDevice );
  
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    GADGET_DEBUG2("Unable to upload samples to GPU memory: %s",
		  cudaGetErrorString(err));
    return GADGET_FAIL;
  }

//GADGET_DEBUG1("Samples uploaded to GPU\n");
  
  return GADGET_OK;
}

int GPUCGGadget::allocate_csm_buffer()
{
  
  if (csm_buffer_dev_ptr_) cudaFree(csm_buffer_dev_ptr_);
  if (csm_acc_coil_image_os_dev_ptr_) cudaFree(csm_acc_coil_image_os_dev_ptr_);

  cudaMalloc( (void**) &csm_buffer_dev_ptr_, 
  	      channels_*prod(matrix_size_os_)*csm_buffer_length_*sizeof(cuFloatComplex) );
  
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    GADGET_DEBUG2("Unable to allocate CSM buffer: %s \n",
		  cudaGetErrorString(err));
    
    return GADGET_FAIL;
  }
  
  clear_image( prod(matrix_size_os_)*channels_*csm_buffer_length_, 
	       make_cuFloatComplex(0.0f, 0.0f), 
	       csm_buffer_dev_ptr_ );

  //Allocate memory for accumulated coil images
  cudaMalloc( (void**) &csm_acc_coil_image_os_dev_ptr_, 
	      prod(matrix_size_os_)*channels_*sizeof(cuFloatComplex) );
  
  err = cudaGetLastError();
  if( err != cudaSuccess ){
    GADGET_DEBUG2("Unable to allocate eccumulated CSM images: %s \n",
		  cudaGetErrorString(err));
    
    return GADGET_FAIL;
  }

  clear_image( prod(matrix_size_os_)*channels_, 
	       make_cuFloatComplex(0.0f, 0.0f), 
	       csm_acc_coil_image_os_dev_ptr_ );

  kspace_acc_coil_images_os_size_ = prod(matrix_size_os_)*channels_;

  //m_csm_needs_reset = true;
  return GADGET_OK;
}
