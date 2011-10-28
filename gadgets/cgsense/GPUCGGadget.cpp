#include "GPUCGGadget.h"
#include "Gadgetron.h"
#include "GadgetMRIHeaders.h"
#include "GadgetXml.h"
#include "ndarray_vector_td_utilities.h"
#include "b1_map.h"

#include "hoNDArray_fileio.h"

#include "tinyxml.h"

GPUCGGadget::GPUCGGadget()
  : slice_no_(0)
  , profiles_per_frame_(32)
  , shared_profiles_(16)
  , channels_(0)
  , samples_per_profile_(0)
  , device_number_(0)
  , number_of_iterations_(5)
  , cg_limit_(1e-6)
  , oversampling_(1.25)
  , kernel_width_(5.5)
  , kappa_(0.1)
  , current_profile_offset_(0)
  , allocated_samples_(0)
  , data_host_ptr_(0x0)
  , is_configured_(false)
  , dcw_computed_(false)
{
  matrix_size_    = uintd2(0,0);
  matrix_size_os_ = uintd2(0,0);
  pass_on_undesired_data_ = true; // We will make one of these for each slice and so data should be passed on.
}

GPUCGGadget::~GPUCGGadget() {}

int GPUCGGadget::process_config( ACE_Message_Block* mb )
{
  GADGET_DEBUG1("GPUCGGadget::process_config\n");

  slice_no_ = get_int_value(std::string("sliceno").c_str());	
  device_number_ = get_int_value(std::string("deviceno").c_str());
  profiles_per_frame_ = get_int_value(std::string("profiles_per_frame").c_str());
  shared_profiles_ = get_int_value(std::string("shared_profiles").c_str());
  number_of_iterations_ = get_int_value(std::string("number_of_iterations").c_str());
  cg_limit_ = get_double_value(std::string("cg_limit").c_str());
  oversampling_ = get_double_value(std::string("oversampling").c_str());
  kernel_width_ = get_double_value(std::string("kernel_width").c_str());
  kappa_ = get_double_value(std::string("kappa").c_str());
  pass_on_undesired_data_ = get_bool_value(std::string("pass_on_undesired_data").c_str());

  if( shared_profiles_ > (profiles_per_frame_>>1) ){
    GADGET_DEBUG1("WARNING: GPUCGGadget::process_config: shared_profiles exceeds half the new samples. Setting to half.\n");
    shared_profiles_ = profiles_per_frame_>>1;
  }

  TiXmlDocument doc;
  doc.Parse(mb->rd_ptr());

  GadgetXMLNode n(&doc);

  if (!is_configured_) {

    cudaDeviceProp deviceProp; 
    if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
      GADGET_DEBUG1( "Error: unable to query device properties.\n" );
      return GADGET_FAIL;
    }

    unsigned int warp_size = deviceProp.warpSize;

    samples_per_profile_ = n.get<long>(std::string("gadgetron.encoding.kspace.readout_length.value"))[0];
    channels_ = n.get<long>(std::string("gadgetron.encoding.channels.value"))[0];


    std::vector<long> dims = n.get<long>(std::string("gadgetron.encoding.kspace.matrix_size.value"));
    matrix_size_ = uintd2(dims[0], dims[1]);

    GADGET_DEBUG2("Matrix size  : [%d,%d] \n", matrix_size_.vec[0], matrix_size_.vec[1]);

    matrix_size_os_ = 
      uintd2(static_cast<unsigned int>(ceil((matrix_size_.vec[0]*oversampling_)/warp_size)*warp_size),
	     static_cast<unsigned int>(ceil((matrix_size_.vec[1]*oversampling_)/warp_size)*warp_size));

    GADGET_DEBUG2("Matrix size OS: [%d,%d] \n", matrix_size_os_.vec[0], matrix_size_os_.vec[1]);

    // Allocate encoding operator for non-Cartesian Sense
    E_ = boost::shared_ptr< cuNonCartesianSenseOperator<float,2> >( new cuNonCartesianSenseOperator<float,2>() );  

    // Allocate preconditioner
    D_ = boost::shared_ptr< cuCGPrecondWeights<float_complext::Type> >( new cuCGPrecondWeights<float_complext::Type>() );

    // Allocate regularization image operator
    R_ = boost::shared_ptr< cuImageOperator<float,float_complext::Type> >( new cuImageOperator<float,float_complext::Type>() );  
    R_->set_weight( kappa_ );

    // Setup solver
    cg_.add_matrix_operator( E_ );  // encoding matrix
    cg_.add_matrix_operator( R_ );  // regularization matrix
    cg_.set_preconditioner ( D_ );  // preconditioning matrix
    cg_.set_iterations( number_of_iterations_ );
    cg_.set_limit( cg_limit_ ); 
    cg_.set_output_mode( cuCGSolver<float, float_complext::Type>::OUTPUT_VERBOSE ); // TODO: once it is all working, change to silent output

    if( configure_channels() == GADGET_FAIL )
      return GADGET_FAIL;
    
    is_configured_ = true;
  }

  return GADGET_OK;
}


int GPUCGGadget::configure_channels()
{
  // We do not have a csm yet, so initialize a dummy one to purely ones
  boost::shared_ptr< cuNDArray<float_complext::Type> > csm = boost::shared_ptr< cuNDArray<float_complext::Type> >( new cuNDArray<float_complext::Type> );
  std::vector<unsigned int> csm_dims = uintd_to_vector<2>(matrix_size_); csm_dims.push_back( channels_ );
  
  if( csm->create( &csm_dims ) == 0x0 ) {
    GADGET_DEBUG1( "Error: unable to create csm.\n" );
    return GADGET_FAIL;
  }
  
  if( !cuNDA_clear<float_complext::Type>( csm.get(), get_one<float_complext::Type>() ) ){
    GADGET_DEBUG1( "Error: unable to clear csm.\n" );
    return GADGET_FAIL;
  }
  
  // Setup matrix operator
  E_->set_csm(csm);
  
  if( E_->setup( matrix_size_, matrix_size_os_, kernel_width_ ) < 0 ){
    GADGET_DEBUG1( "Error: unable to setup encoding operator.\n" );
    return GADGET_FAIL;
  }
  
  // Allocate rhs buffer
  rhs_buffer_ = boost::shared_ptr< cuSenseRHSBuffer<float,2> >( new cuSenseRHSBuffer<float,2>() );
  rhs_buffer_->set_num_coils( channels_ );
  rhs_buffer_->set_sense_operator( E_ );

  return GADGET_OK;
}

int GPUCGGadget::process(GadgetContainerMessage<GadgetMessageAcquisition>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
  if (!is_configured_) {
    GADGET_DEBUG1("Data received before configuration complete\n");
    return GADGET_FAIL;
  }

  //Is this data for me?
  if (m1->getObjectPtr()->idx.slice != slice_no_) {

    //This data is not for me
    if (pass_on_undesired_data_) {
      this->next()->putq(m1);
    } else {
      GADGET_DEBUG2("Dropping slice: %d\n", m1->getObjectPtr()->idx.slice);
      m1->release();
    }
    return GADGET_OK;
  }

  // Check if some upstream gadget has modified the number of channels or samples per profile 
  // since the global configuration is no longer valid then...
  //

  if( m1->getObjectPtr()->samples != samples_per_profile_ ) {
    GADGET_DEBUG2("Adjusting #samples per profile from %d to %d", samples_per_profile_,  m1->getObjectPtr()->samples );
    samples_per_profile_ = m1->getObjectPtr()->samples;
    allocated_samples_ = 0; // the samples buffers are freed and re-allocated in 'upload_samples()'
  }

  if( m1->getObjectPtr()->channels != channels_ ) {
    GADGET_DEBUG2("Adjusting #channels from %d to %d", channels_,  m1->getObjectPtr()->channels );
    channels_ = m1->getObjectPtr()->channels;
    allocated_samples_ = 0; // the samples buffers are freed and re-allocated in 'upload_samples()'
    if( configure_channels() == GADGET_FAIL ) // Update buffers dependant on #channels
      return GADGET_FAIL;    
  }

  buffer_.enqueue_tail(m1);

  if ((int)buffer_.message_count() >= profiles_per_frame_) {

    boost::shared_ptr< cuNDArray<floatd2::Type> > traj = calculate_trajectory();

    if ( traj.get() == 0x0 ) {
      GADGET_DEBUG1("Failed to calculate trajectory\n");
      return GADGET_FAIL;
    }

    boost::shared_ptr< cuNDArray<float> > dcw;
    if( !dcw_computed_){
      dcw = calculate_density_compensation();
      if( dcw.get() == 0x0 ) {
	GADGET_DEBUG1("Failed to calculate density compensation\n");
	return GADGET_FAIL;
      }
      E_->set_dcw(dcw);
      dcw_computed_ = true;
    }

    boost::shared_ptr< cuNDArray<float_complext::Type> > device_samples = upload_samples();
    if( device_samples == 0x0 ) {
      GADGET_DEBUG1("Failed to upload samples to the GPU\n");
      return GADGET_FAIL;
    }

    if( E_->preprocess(traj.get()) < 0 ) {
      GADGET_DEBUG1("Error during cgOperatorNonCartesianSense::preprocess()\n");
      return GADGET_FAIL;
    }

    rhs_buffer_->add_frame_data( device_samples.get(), traj.get() );

    boost::shared_ptr< cuNDArray<float_complext::Type> > csm_data = rhs_buffer_->get_acc_coil_images();
    if( !csm_data.get() ){
      GADGET_DEBUG1("Error during accumulation buffer computation\n");
      return GADGET_FAIL;
    }

    // Estimate CSM
    boost::shared_ptr< cuNDArray<float_complext::Type> > csm = estimate_b1_map<float,2>( csm_data.get() );
    E_->set_csm(csm);

    boost::shared_ptr< std::vector<unsigned int> > reg_dims = csm_data->get_dimensions();
    reg_dims->pop_back();

    cuNDArray<float_complext::Type> reg_image;
    if( reg_image.create(reg_dims.get()) == 0x0 ){
      GADGET_DEBUG1("Error allocating regularization image on device\n");
      return GADGET_FAIL;
    }

    if( E_->mult_csm_conj_sum( csm_data.get(), &reg_image ) < 0 ){
      GADGET_DEBUG1("Error combining coils to regularization image\n");
      return GADGET_FAIL;
    }
	
    R_->compute(&reg_image);

    // TODO: error check these computations

    // Define preconditioning weights
    boost::shared_ptr< cuNDArray<float> > _precon_weights = cuNDA_ss<float,float_complext::Type>( csm.get(), 2 );
    cuNDA_axpy<float>( kappa_, R_->get(), _precon_weights.get() );  
    cuNDA_reciprocal_sqrt<float>( _precon_weights.get() );
    boost::shared_ptr< cuNDArray<float_complext::Type> > precon_weights = cuNDA_real_to_complext<float>( _precon_weights.get() );
    _precon_weights.reset();
    D_->set_weights( precon_weights );

    // Form rhs
    std::vector<unsigned int> rhs_dims = uintd_to_vector<2>(matrix_size_);
    cuNDArray<float_complext::Type> rhs; 
		
    if( rhs.create(&rhs_dims) == 0x0 ){
      GADGET_DEBUG1("failed to create rhs\n");
      return GADGET_FAIL;
    }

    if( E_->mult_MH( device_samples.get(), &rhs ) < 0 ){
      GADGET_DEBUG1("failed to compute rhs\n");
      return GADGET_FAIL;
    }

    boost::shared_ptr< cuNDArray<float_complext::Type> > cgresult = cg_.solve(&rhs);

    if (!cgresult.get()) {
      GADGET_DEBUG1("iterative_sense_compute failed\n");
      return GADGET_FAIL;
    }
	
    //Now pass the reconstructed image on
    GadgetContainerMessage<GadgetMessageImage>* cm1 = 
      new GadgetContainerMessage<GadgetMessageImage>();

    GadgetContainerMessage< hoNDArray< std::complex<float> > >* cm2 = 
      new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

    cm1->cont(cm2);

    std::vector<unsigned int> img_dims(2);
    img_dims[0] = matrix_size_.vec[0];
    img_dims[1] = matrix_size_.vec[1];

    if (cm2->getObjectPtr()->create(&img_dims) == 0x0) {
      GADGET_DEBUG1("Unable to allocate host image array");
      cm1->release();
      return GADGET_FAIL;
    }

    size_t data_length = prod(matrix_size_);

    cudaMemcpy(cm2->getObjectPtr()->get_data_ptr(),
	       cgresult->get_data_ptr(),
	       data_length*sizeof(std::complex<float>),
	       cudaMemcpyDeviceToHost);			

    cudaError_t err = cudaGetLastError();
    if( err != cudaSuccess ){
      GADGET_DEBUG2("Unable to copy result from device to host: %s", cudaGetErrorString(err));
      cm1->release();
      return GADGET_FAIL;
    }

    cm1->getObjectPtr()->matrix_size[0] = img_dims[0];
    cm1->getObjectPtr()->matrix_size[1] = img_dims[1];
    cm1->getObjectPtr()->matrix_size[2] = 1;
    cm1->getObjectPtr()->channels       = 1;
    cm1->getObjectPtr()->data_idx_min       = m1->getObjectPtr()->min_idx;
    cm1->getObjectPtr()->data_idx_max       = m1->getObjectPtr()->max_idx;
    cm1->getObjectPtr()->data_idx_current   = m1->getObjectPtr()->idx;	

    memcpy(cm1->getObjectPtr()->position,m1->getObjectPtr()->position, sizeof(float)*3);
    memcpy(cm1->getObjectPtr()->quarternion,m1->getObjectPtr()->quarternion, sizeof(float)*4);

    if (this->next()->putq(cm1) < 0) {
      GADGET_DEBUG1("Failed to result image on to Q\n");
      cm1->release();
      return GADGET_FAIL;
    }

    //Dequeue the message we don't need anymore
    ACE_Message_Block* mb_tmp;
    for (int i = 0; i < (profiles_per_frame_-shared_profiles_); i++) {
      buffer_.dequeue_head(mb_tmp);
      mb_tmp->release();
      current_profile_offset_++; 
    }
  }

  return GADGET_OK;
}


int GPUCGGadget::copy_samples_for_profile(float* host_base_ptr,
					  std::complex<float>* data_base_ptr,
					  int profile_no,
					  int channel_no)
{

  memcpy(host_base_ptr + 
	 (channel_no*allocated_samples_ + profile_no*samples_per_profile_) * 2,
	 data_base_ptr + channel_no*samples_per_profile_, 
	 sizeof(float)*samples_per_profile_*2);

  return GADGET_OK;
}

boost::shared_ptr< cuNDArray<float_complext::Type> >  GPUCGGadget::upload_samples()
{

  int samples_needed = 
    samples_per_profile_*
    profiles_per_frame_;

  if (samples_needed != allocated_samples_) {
    
    if( data_host_ptr_ ){
      delete[] data_host_ptr_;
      data_host_ptr_ = 0x0;
      allocated_samples_ = 0;
    }

    try {
      data_host_ptr_ = new float[channels_*samples_needed*2];
    } catch (...) {
      GADGET_DEBUG1("Failed to allocate host memory for samples\n");
      return boost::shared_ptr< cuNDArray<float_complext::Type> >();
    }

    allocated_samples_ = samples_needed;

  }

  ACE_Message_Queue_Reverse_Iterator<ACE_MT_SYNCH> it(buffer_);
  int profiles_copied = 0;
  GadgetContainerMessage<GadgetMessageAcquisition>* m1;
  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2;
  ACE_Message_Block* mb;

  while (profiles_copied < profiles_per_frame_) {
    it.next(mb);

    m1 = dynamic_cast< GadgetContainerMessage< GadgetMessageAcquisition >* >(mb);
    if (!m1) {
      GADGET_DEBUG1("Failed to dynamic cast message\n");
      return boost::shared_ptr< cuNDArray<float_complext::Type> >();
    }

    m2 = dynamic_cast< GadgetContainerMessage< hoNDArray< std::complex<float> > >* > (m1->cont());

    if (!m2) {
      GADGET_DEBUG1("Failed to dynamic cast message\n");
      return boost::shared_ptr< cuNDArray<float_complext::Type> >();
    }

    std::complex<float> *d = m2->getObjectPtr()->get_data_ptr();
    int current_profile = profiles_per_frame_-profiles_copied-1;

    for (int i = 0; i < channels_; i++) {
      copy_samples_for_profile( data_host_ptr_, d, current_profile, i );
    }

    it.advance();   
    profiles_copied++;
  }
	
  std::vector<unsigned int> dims; dims.push_back(samples_needed); dims.push_back(channels_);
  hoNDArray<float_complext::Type> tmp;
  if( tmp.create( &dims, (float_complext::Type*)data_host_ptr_, false ) == 0x0 ){
    GADGET_DEBUG1("Failed to create temporary host data array\n");
    return boost::shared_ptr< cuNDArray<float_complext::Type> >();
  }

  boost::shared_ptr< cuNDArray<float_complext::Type> > device_samples ( new cuNDArray<float_complext::Type>(&tmp) );
	
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    GADGET_DEBUG2("Unable to upload samples to GPU memory: %s", cudaGetErrorString(err));
    return boost::shared_ptr< cuNDArray<float_complext::Type> >();
  }
  
  return device_samples;
}
