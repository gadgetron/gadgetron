#include "gpuSpiralSensePrepGadget.h"
#include "SenseJob.h"
#include "Gadgetron.h"
#include "cuNDArray_utils.h"
#include "vector_td_utilities.h"
#include "hoNDArray_fileio.h"
#include "vector_td.h"
#include "vector_td_operators.h"
#include "check_CUDA.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "vds.h"

#include <algorithm>
#include <vector>

namespace Gadgetron{

  gpuSpiralSensePrepGadget::gpuSpiralSensePrepGadget()
    : samples_to_skip_start_(0)
    , samples_to_skip_end_(0)
    , samples_per_interleave_(0)
    , prepared_(false)
    , use_multiframe_grouping_(false)
    , acceleration_factor_(0)
  {
    GADGET_DEBUG1("Initializing Spiral\n");
    set_parameter(std::string("buffer_using_solver").c_str(), "false");
    set_parameter(std::string("propagate_csm_from_set").c_str(), "-1");
    set_parameter(std::string("buffer_convolution_kernel_width").c_str(), "5.5");
    set_parameter(std::string("buffer_convolution_oversampling_factor").c_str(), "1.25");
    set_parameter(std::string("reconstruction_os_factor_x").c_str(), "1.0");
    set_parameter(std::string("reconstruction_os_factor_y").c_str(), "1.0");
  }

  gpuSpiralSensePrepGadget::~gpuSpiralSensePrepGadget() {}

  int gpuSpiralSensePrepGadget::process_config(ACE_Message_Block* mb)
  {

    int number_of_devices = 0;
    if (cudaGetDeviceCount(&number_of_devices)!= cudaSuccess) {
      GADGET_DEBUG1( "Error: unable to query number of CUDA devices.\n" );
      return GADGET_FAIL;
    }

    if (number_of_devices == 0) {
      GADGET_DEBUG1( "Error: No available CUDA devices.\n" );
      return GADGET_FAIL;
    }

    device_number_ = get_int_value(std::string("deviceno").c_str());

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

    propagate_csm_from_set_ = get_int_value(std::string("propagate_csm_from_set").c_str());

    if( propagate_csm_from_set_ > 0 ){
      GADGET_DEBUG2("Currently, only set 0 can propagate coil sensitivity maps. Set %d was specified.\n", propagate_csm_from_set_ );
      return GADGET_FAIL;
    }

    if( propagate_csm_from_set_ >= 0 ){
      GADGET_DEBUG2("Propagating csm from set %d to all sets\n", propagate_csm_from_set_ );
    }

    buffer_using_solver_ = get_bool_value(std::string("buffer_using_solver").c_str());
    use_multiframe_grouping_ = get_bool_value(std::string("use_multiframe_grouping").c_str());

    if( buffer_using_solver_ && !use_multiframe_grouping_ ){
      GADGET_DEBUG1("Enabling 'buffer_using_solver' requires also enabling 'use_multiframe_grouping'.\n" );
      return GADGET_FAIL;
    }

    // Start parsing the ISMRMRD XML header
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

    //ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    // Determine reconstruction matrix sizes
    //

    kernel_width_ = get_double_value(std::string("buffer_convolution_kernel_width").c_str());
    oversampling_factor_ = get_double_value(std::string("buffer_convolution_oversampling_factor").c_str());
    
    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(r_space.matrixSize().x()*get_double_value(std::string("reconstruction_os_factor_x").c_str())))+warp_size-1)/warp_size)*warp_size);  
    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(r_space.matrixSize().y()*get_double_value(std::string("reconstruction_os_factor_y").c_str())))+warp_size-1)/warp_size)*warp_size);
      
    image_dimensions_recon_os_ = uintd2
      (((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
       ((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);
    
    // In case the warp_size constraint kicked in
    oversampling_factor_ = float(image_dimensions_recon_os_[0])/float(image_dimensions_recon_[0]);
    
    //
    // Setup the spiral trajectory
    //

    if (!(*e_seq.begin()).trajectoryDescription().present()) {
      GADGET_DEBUG1("Trajectory description needed to calculate trajectory");
      return GADGET_FAIL;
    }

    ISMRMRD::trajectoryDescriptionType traj_desc = (*e_seq.begin()).trajectoryDescription().get();

    if (std::strcmp(traj_desc.identifier().c_str(), "HargreavesVDS2000")) {
      GADGET_DEBUG1("Expected trajectory description identifier 'HargreavesVDS2000', not found.");
      return GADGET_FAIL;
    }

    long interleaves = -1;
    long fov_coefficients = -1;
    long sampling_time_ns = -1;
    double max_grad = -1.0;
    double max_slew = -1.0;
    double fov_coeff = -1.0;
    double kr_max = -1.0;

    for (ISMRMRD::trajectoryDescriptionType::userParameterLong_sequence::iterator i (traj_desc.userParameterLong().begin ()); i != traj_desc.userParameterLong().end(); ++i) {
      if (std::strcmp(i->name().c_str(),"interleaves") == 0) {
	interleaves = i->value();
      } else if (std::strcmp(i->name().c_str(),"fov_coefficients") == 0) {
	fov_coefficients = i->value();
      } else if (std::strcmp(i->name().c_str(),"SamplingTime_ns") == 0) {
	sampling_time_ns = i->value();
      } else {
	GADGET_DEBUG2("WARNING: unused trajectory parameter %s found\n", i->name().c_str());
      }
    }

    for (ISMRMRD::trajectoryDescriptionType::userParameterDouble_sequence::iterator i (traj_desc.userParameterDouble().begin ()); i != traj_desc.userParameterDouble().end(); ++i) {
      if (std::strcmp(i->name().c_str(),"MaxGradient_G_per_cm") == 0) {
	max_grad = i->value();
      } else if (std::strcmp(i->name().c_str(),"MaxSlewRate_G_per_cm_per_s") == 0) {
	max_slew = i->value();
      } else if (std::strcmp(i->name().c_str(),"FOVCoeff_1_cm") == 0) {
	fov_coeff = i->value();
      } else if (std::strcmp(i->name().c_str(),"krmax_per_cm") == 0) {
	kr_max= i->value();
      } else {
	GADGET_DEBUG2("WARNING: unused trajectory parameter %s found\n", i->name().c_str());
      }
    }

    if ((interleaves < 0) || (fov_coefficients < 0) || (sampling_time_ns < 0) || (max_grad < 0) || (max_slew < 0) || (fov_coeff < 0) || (kr_max < 0)) {
      GADGET_DEBUG1("Appropriate parameters for calculating spiral trajectory not found in XML configuration\n");
      return GADGET_FAIL;
    }

    Tsamp_ns_ = sampling_time_ns;
    Nints_ = interleaves;
    interleaves_ = static_cast<int>(Nints_);

    gmax_ = max_grad;
    smax_ = max_slew;
    krmax_ = kr_max;
    fov_ = fov_coeff;

    samples_to_skip_start_  = 0; //n.get<int>(std::string("samplestoskipstart.value"))[0];
    samples_to_skip_end_    = -1; //n.get<int>(std::string("samplestoskipend.value"))[0];

    fov_vec_.push_back(r_space.fieldOfView_mm().x());
    fov_vec_.push_back(r_space.fieldOfView_mm().y());
    fov_vec_.push_back(r_space.fieldOfView_mm().z());

    slices_ = e_limits.slice().present() ? e_limits.slice().get().maximum() + 1 : 1;
    sets_ = e_limits.set().present() ? e_limits.set().get().maximum() + 1 : 1;

    buffer_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);

    image_headers_queue_ = 
      boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);

    size_t bsize = sizeof(GadgetContainerMessage<ISMRMRD::ImageHeader>)*100*Nints_;

    for( unsigned int i=0; i<slices_*sets_; i++ ){
      image_headers_queue_[i].high_water_mark(bsize);
      image_headers_queue_[i].low_water_mark(bsize);
    }

    GADGET_DEBUG2("smax:                    %f\n", smax_);
    GADGET_DEBUG2("gmax:                    %f\n", gmax_);
    GADGET_DEBUG2("Tsamp_ns:                %d\n", Tsamp_ns_);
    GADGET_DEBUG2("Nints:                   %d\n", Nints_);
    GADGET_DEBUG2("fov:                     %f\n", fov_);
    GADGET_DEBUG2("krmax:                   %f\n", krmax_);
    GADGET_DEBUG2("samples_to_skip_start_ : %d\n", samples_to_skip_start_);
    GADGET_DEBUG2("samples_to_skip_end_   : %d\n", samples_to_skip_end_);
    GADGET_DEBUG2("recon matrix_size_x    : %d\n", image_dimensions_recon_[0]);
    GADGET_DEBUG2("recon matrix_size_y    : %d\n", image_dimensions_recon_[1]);

    return GADGET_OK;
  }

  int gpuSpiralSensePrepGadget::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2)
  {
    // Noise should have been consumed by the noise adjust, but just in case...
    //

    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    if (is_noise) {
      m1->release();
      return GADGET_OK;
    }

    if (!prepared_) {

      int     nfov   = 1;         /*  number of fov coefficients.             */
      int     ngmax  = 1e5;       /*  maximum number of gradient samples      */
      double  *xgrad;             /*  x-component of gradient.                */
      double  *ygrad;             /*  y-component of gradient.                */
      double  *x_trajectory;
      double  *y_trajectory;
      double  *weighting;
      int     ngrad;
      //int     count;
      double sample_time = (1.0*Tsamp_ns_) * 1e-9;

      /*	call c-function here to calculate gradients */
      calc_vds(smax_,gmax_,sample_time,sample_time,Nints_,&fov_,nfov,krmax_,ngmax,&xgrad,&ygrad,&ngrad);
      samples_per_interleave_ = std::min(ngrad,static_cast<int>(m1->getObjectPtr()->number_of_samples));

      GADGET_DEBUG2("Using %d samples per interleave\n", samples_per_interleave_);

      /* Calcualte the trajectory and weights*/
      calc_traj(xgrad, ygrad, samples_per_interleave_, Nints_, sample_time, krmax_, &x_trajectory, &y_trajectory, &weighting);

      host_traj_ = boost::shared_ptr< hoNDArray<floatd2> >(new hoNDArray<floatd2>);
      host_weights_ = boost::shared_ptr< hoNDArray<float> >(new hoNDArray<float>);

      std::vector<unsigned int> trajectory_dimensions;
      trajectory_dimensions.push_back(samples_per_interleave_*Nints_);

      host_traj_->create(&trajectory_dimensions);
      host_weights_->create(&trajectory_dimensions);

      {
	float* co_ptr = reinterpret_cast<float*>(host_traj_->get_data_ptr());
	float* we_ptr =  reinterpret_cast<float*>(host_weights_->get_data_ptr());
	
	for (int i = 0; i < (samples_per_interleave_*Nints_); i++) {
	  co_ptr[i*2]   = -x_trajectory[i]/2;
	  co_ptr[i*2+1] = -y_trajectory[i]/2;
	  we_ptr[i] = weighting[i];
	}
      }

      delete [] xgrad;
      delete [] ygrad;
      delete [] x_trajectory;
      delete [] y_trajectory;
      delete [] weighting;

      // Setup the NFFT plan
      //

      cuNDArray<floatd2> traj(*host_traj_);
      dcw_buffer_ = boost::shared_ptr< cuNDArray<float> >( new cuNDArray<float>(*host_weights_) );
	
      nfft_plan_.setup( from_std_vector<unsigned int,2>(image_dimensions_recon_), image_dimensions_recon_os_, kernel_width_ );
      nfft_plan_.preprocess(&traj, cuNFFT_plan<float,2>::NFFT_PREP_NC2C);

      // Setup the non-Cartesian Sense encoding operator 
      //
      
      E_ = boost::shared_ptr< cuNonCartesianSenseOperator<float,2> >(new cuNonCartesianSenseOperator<float,2>);
      E_->setup( from_std_vector<unsigned int,2>(image_dimensions_recon_), image_dimensions_recon_os_, kernel_width_ );
      
      // Setup cg solver if the csm/regularization image is to be based hereon
      //

      if( buffer_using_solver_ ){

	E_->set_dcw(dcw_buffer_);

	D_ = boost::shared_ptr< cuCgPreconditioner<float_complext> >( new cuCgPreconditioner<float_complext>() );
	cg_.set_encoding_operator( E_ );
	cg_.set_preconditioner( D_ );
	cg_.set_max_iterations( 2 );
	cg_.set_tc_tolerance( 1e-6 );
	cg_.set_output_mode( cuCgSolver<float_complext>::OUTPUT_SILENT);
      }

      prepared_ = true;
    }

    // Allocate host data buffer if it is NULL
    //

    if (!host_data_buffer_.get()) {

      std::vector<unsigned int> data_dimensions;
      data_dimensions.push_back(samples_per_interleave_*interleaves_);
      data_dimensions.push_back(m1->getObjectPtr()->active_channels);

      host_data_buffer_ = boost::shared_array< hoNDArray<float_complext> >
	(new hoNDArray<float_complext>[slices_*sets_]);
      
      if (!host_data_buffer_.get()) {
	GADGET_DEBUG1("Unable to allocate array for host data buffer\n");
	return GADGET_FAIL;
      }

      for (unsigned int i = 0; i < slices_*sets_; i++) {
	host_data_buffer_[i].create(&data_dimensions);
	host_data_buffer_[i].fill(0.0f);
      }
    }

    // Allocate various counters if they are NULL
    //

    if( !image_counter_.get() ){
      image_counter_ = boost::shared_array<long>(new long[slices_*sets_]);
      for( unsigned int i=0; i<slices_*sets_; i++ )
	image_counter_[i] = 0;
    }

    if( !interleaves_counter_singleframe_.get() ){
      interleaves_counter_singleframe_ = boost::shared_array<long>(new long[slices_*sets_]);
      for( unsigned int i=0; i<slices_*sets_; i++ )
	interleaves_counter_singleframe_[i] = 0;
    }

    if( !interleaves_counter_multiframe_.get() ){
      interleaves_counter_multiframe_ = boost::shared_array<long>(new long[slices_*sets_]);
      for( unsigned int i=0; i<slices_*sets_; i++ )
	interleaves_counter_multiframe_[i] = 0;
    }

    // Define some utility variables
    //

    unsigned int samples_to_copy = m1->getObjectPtr()->number_of_samples-samples_to_skip_end_;
    unsigned int interleave = m1->getObjectPtr()->idx.kspace_encode_step_1;
    unsigned int slice = m1->getObjectPtr()->idx.slice;
    unsigned int set = m1->getObjectPtr()->idx.set;
    unsigned int samples_per_channel =  host_data_buffer_[set*slices_+slice].get_size(0);

    // Some book-keeping to keep track of the frame count
    //

    interleaves_counter_singleframe_[set*slices_+slice]++;
    interleaves_counter_multiframe_[set*slices_+slice]++;

    // Duplicate the profile to avoid double deletion in case problems are encountered below.
    // Enque profile until all profiles for the reconstruction have been received.
    //
    
    buffer_[set*slices_+slice].enqueue_tail(duplicate_profile(m1));
    
    // Copy profile into the accumulation buffer for csm/regularization estimation
    //

    ISMRMRD::AcquisitionHeader base_head = *m1->getObjectPtr();

    if (samples_to_skip_end_ == -1) {
      samples_to_skip_end_ = m1->getObjectPtr()->number_of_samples-samples_per_interleave_;
      GADGET_DEBUG2("Adjusting samples_to_skip_end_ = %d\n", samples_to_skip_end_);
    }

    std::complex<float>* data_ptr = reinterpret_cast< std::complex<float>* >
      (host_data_buffer_[set*slices_+slice].get_data_ptr());

    std::complex<float>* profile_ptr = m2->getObjectPtr()->get_data_ptr();

    for (unsigned int c = 0; c < m1->getObjectPtr()->active_channels; c++) {
      memcpy(data_ptr+c*samples_per_channel+interleave*samples_to_copy,
	     profile_ptr+c*m1->getObjectPtr()->number_of_samples, samples_to_copy*sizeof(std::complex<float>));
    }

    // Have we received sufficient data for a new frame?
    //

    bool is_last_scan_in_slice = 
      ISMRMRD::FlagBit(ISMRMRD::ACQ_LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags);

    if (is_last_scan_in_slice) {

      // This was the final profile of a frame
      //

      if( Nints_%interleaves_counter_singleframe_[set*slices_+slice] ){
	GADGET_DEBUG1("Unexpected number of interleaves encountered in frame\n");
	return GADGET_FAIL;
      }

      // Has the acceleration factor changed?
      //

      if( acceleration_factor_ != Nints_/interleaves_counter_singleframe_[set*slices_+slice] ){

	GADGET_DEBUG1("Change of acceleration factor detected\n");
	acceleration_factor_ =  Nints_/interleaves_counter_singleframe_[set*slices_+slice];

	// The encoding operator needs to have its domain/codomain dimensions set accordingly
	//
	
	if( buffer_using_solver_ ){

	  std::vector<unsigned int> domain_dims = image_dimensions_recon_;
	  
	  std::vector<unsigned int> codomain_dims = *host_traj_->get_dimensions();
	  codomain_dims.push_back(m1->getObjectPtr()->active_channels);
	  
	  E_->set_domain_dimensions(&domain_dims);
	  E_->set_codomain_dimensions(&codomain_dims);

	  cuNDArray<floatd2> traj(*host_traj_);
	  E_->preprocess(&traj);
	}
      }

      // Prepare an image header for this frame
      //

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
      header->getObjectPtr()->matrix_size[2] = acceleration_factor_;

      header->getObjectPtr()->field_of_view[0] = fov_vec_[0];
      header->getObjectPtr()->field_of_view[1] = fov_vec_[1];
      header->getObjectPtr()->field_of_view[2] = fov_vec_[2];

      header->getObjectPtr()->channels = base_head->active_channels;
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
      header->getObjectPtr()->image_index = image_counter_[set*slices_+slice]++; 
      header->getObjectPtr()->image_series_index = set*slices_+slice;

      // Enque header until we are ready to assemble a Sense job
      //

      image_headers_queue_[set*slices_+slice].enqueue_tail(header);

      // Check if it is time to reconstruct.
      // I.e. prepare and pass a Sense job downstream...
      //

      if( !use_multiframe_grouping_ || 
	  (use_multiframe_grouping_ && interleaves_counter_multiframe_[set*slices_+slice] == Nints_) ){

	unsigned int num_coils = m1->getObjectPtr()->active_channels;
	
	// Compute coil images from the fully sampled data buffer
	//

	std::vector<unsigned int> image_dims;
	image_dims.push_back(image_dimensions_recon_[0]);
	image_dims.push_back(image_dimensions_recon_[1]);
	image_dims.push_back(num_coils);
	
	cuNDArray<float_complext> image(&image_dims);
	cuNDArray<float_complext> data(&host_data_buffer_[set*slices_+slice]);
	
	nfft_plan_.compute( &data, &image, dcw_buffer_.get(), cuNFFT_plan<float,2>::NFFT_BACKWARDS_NC2C );

	// Check if we need to compute a new csm
	//
	
	if( propagate_csm_from_set_ < 0 || propagate_csm_from_set_ == set ){	  	  
	  csm_ = estimate_b1_map<float,2>( &image ); // Estimates csm
	}
	else{
	  //GADGET_DEBUG2("Set %d is reusing the csm from set %d\n", set, propagate_csm_from_set_);
	  if( csm_.get() == 0x0 ){
	    GADGET_DEBUG1("Error, csm has not been computed\n");
	    return GADGET_FAIL;
	  }	  
	}
	E_->set_csm(csm_);

	// Compute regularization using basic coil combination
	//
	
	image_dims.pop_back();
	cuNDArray<float_complext> reg_image(&image_dims);
	E_->mult_csm_conj_sum( &image, &reg_image );
	
	if( buffer_using_solver_ ){
	  
	  // Compute regularization using cg solver
	  //
	  
	  // Define preconditioning weights
	  boost::shared_ptr< cuNDArray<float> > _precon_weights = sum(abs_square(csm_.get()).get(), 2);
	  reciprocal_sqrt_inplace(_precon_weights.get());	
	  boost::shared_ptr< cuNDArray<float_complext> > precon_weights = real_to_complex<float_complext>( _precon_weights.get() );
	  _precon_weights.reset();
	  D_->set_weights( precon_weights );
	  
	  // Solve from the plain coil combination
	  reg_image = *cg_.solve_from_rhs(&reg_image);
	}

	// Get ready to fill in the Sense job
	//

	boost::shared_ptr< hoNDArray<float_complext> > csm_host = csm_->to_host();
	boost::shared_ptr< hoNDArray<float_complext> > reg_host = reg_image.to_host();

	unsigned int profiles_buffered = buffer_[set*slices_+slice].message_count();

	std::vector<unsigned int> ddimensions;
	ddimensions.push_back(samples_per_interleave_*interleaves_counter_singleframe_[set*slices_+slice]*
			      ((use_multiframe_grouping_) ? acceleration_factor_ : 1));
	ddimensions.push_back(num_coils);
	
	boost::shared_ptr< hoNDArray<float_complext> > data_host(new hoNDArray<float_complext>(&ddimensions));

	ddimensions.clear();
	ddimensions.push_back(samples_per_interleave_*interleaves_counter_singleframe_[set*slices_+slice]);
	ddimensions.push_back((use_multiframe_grouping_) ? acceleration_factor_ : 1);

	boost::shared_ptr< hoNDArray<floatd2> > traj_host(new hoNDArray<floatd2>(&ddimensions));
	boost::shared_ptr< hoNDArray<float> > dcw_host(new hoNDArray<float>(&ddimensions));
	
	for (unsigned int p = 0; p < profiles_buffered; p++) {
	  ACE_Message_Block* mbq;
	  if (buffer_[set*slices_+slice].dequeue_head(mbq) < 0) {
	    GADGET_DEBUG1("Message dequeue failed\n");
	    return GADGET_FAIL;
	  }

	  GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* acq = 
	    AsContainerMessage<ISMRMRD::AcquisitionHeader>(mbq);

	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* daq = 
	    AsContainerMessage<hoNDArray< std::complex<float> > >(mbq->cont());

	  if (!acq || !daq) {
	    GADGET_DEBUG1("Unable to interpret data on message Q\n");
	    return GADGET_FAIL;
	  }

	  for (unsigned int c = 0; c < num_coils; c++) {
	    float_complext* data_ptr = data_host->get_data_ptr();
	    data_ptr += c*samples_per_interleave_*profiles_buffered+p*samples_per_interleave_;

	    std::complex<float>* r_ptr = daq->getObjectPtr()->get_data_ptr();
	    r_ptr += c*daq->getObjectPtr()->get_size(0);

	    memcpy(data_ptr,r_ptr,samples_per_interleave_*sizeof(float_complext));
	  }

	  floatd2* traj_ptr = traj_host->get_data_ptr();
	  traj_ptr += p*samples_per_interleave_;

	  floatd2* t_ptr = host_traj_->get_data_ptr();
	  t_ptr += acq->getObjectPtr()->idx.kspace_encode_step_1*samples_per_interleave_;

	  memcpy(traj_ptr,t_ptr,samples_per_interleave_*sizeof(floatd2));

	  float* dcw_ptr = dcw_host->get_data_ptr();
	  dcw_ptr += p*samples_per_interleave_;

	  float* d_ptr = host_weights_->get_data_ptr();
	  d_ptr += acq->getObjectPtr()->idx.kspace_encode_step_1*samples_per_interleave_;

	  memcpy(dcw_ptr,d_ptr,samples_per_interleave_*sizeof(float));

	  mbq->release();
	}

	GadgetContainerMessage< SenseJob >* m4 = new GadgetContainerMessage< SenseJob >();

	m4->getObjectPtr()->dat_host_ = data_host;
	m4->getObjectPtr()->csm_host_ = csm_host;
	m4->getObjectPtr()->reg_host_ = reg_host;
	m4->getObjectPtr()->tra_host_ = traj_host;
	m4->getObjectPtr()->dcw_host_ = dcw_host;

	// Pull the image headers out of the queue
	//
	
	long frames_per_reconstruction = (use_multiframe_grouping_) ? acceleration_factor_ : 1;
      
	if( image_headers_queue_[set*slices_+slice].message_count() != frames_per_reconstruction ){
	  m4->release();
	  GADGET_DEBUG2("Unexpected size of image header queue: %d, %d\n", 
			image_headers_queue_[set*slices_+slice].message_count(), frames_per_reconstruction);
	  return GADGET_FAIL;
	}
	
	m4->getObjectPtr()->image_headers_ =
	  boost::shared_array<ISMRMRD::ImageHeader>( new ISMRMRD::ImageHeader[frames_per_reconstruction] );
	
	for( unsigned int i=0; i<frames_per_reconstruction; i++ ){	
	  
	  ACE_Message_Block *mbq;
	  
	  if( image_headers_queue_[set*slices_+slice].dequeue_head(mbq) < 0 ) {
	    m4->release();
	    GADGET_DEBUG1("Image header dequeue failed\n");
	    return GADGET_FAIL;
	  }
	  
	  GadgetContainerMessage<ISMRMRD::ImageHeader> *m = AsContainerMessage<ISMRMRD::ImageHeader>(mbq);
	  m4->getObjectPtr()->image_headers_[i] = *m->getObjectPtr();
	  m->release();
	}

	// The Sense Job needs an image header as well. 
	// Let us just copy the initial one...
	
	GadgetContainerMessage<ISMRMRD::ImageHeader> *m3 = new GadgetContainerMessage<ISMRMRD::ImageHeader>;
	*m3->getObjectPtr() = m4->getObjectPtr()->image_headers_[0];
	m3->cont(m4);
	
	if (this->next()->putq(m3) < 0) {
	  GADGET_DEBUG1("Failed to put job on queue.\n");
	  m3->release();
	  return GADGET_FAIL;
	}
	interleaves_counter_multiframe_[set*slices_+slice] = 0;
      }
      interleaves_counter_singleframe_[set*slices_+slice] = 0;
    }
    m1->release();
    return GADGET_OK;
  }

  GadgetContainerMessage<ISMRMRD::AcquisitionHeader>*
  gpuSpiralSensePrepGadget::duplicate_profile( GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *profile )
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

  GADGET_FACTORY_DECLARE(gpuSpiralSensePrepGadget)
}
