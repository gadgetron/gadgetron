#include "gpuSpiralSensePrepGadget.h"
#include "SenseJob.h"
#include "Gadgetron.h"
#include "cuNDArray.h"
#include "hoNDArray_fileio.h"
#include "vector_td.h"
#include "vector_td_operators.h"
#include "check_CUDA.h"
#include "b1_map.h"
#include "cuNonCartesianSenseOperator.h"
#include "GPUTimer.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "vds.h"

#include <algorithm>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace Gadgetron{

  gpuSpiralSensePrepGadget::gpuSpiralSensePrepGadget()
    : samples_to_skip_start_(0)
    , samples_to_skip_end_(0)
    , samples_per_interleave_(0)
    , host_data_buffer_(0)
    , image_counter_(0)
    , image_series_(0)
    , prepared_(false)
    , use_multiframe_grouping_(false)
    , interleaves_counter_singleframe_(0)
    , interleaves_counter_multiframe_(0)
    , acceleration_factor_(0)
  {
    GADGET_DEBUG1("Initializing Spiral\n");
  }

  gpuSpiralSensePrepGadget::~gpuSpiralSensePrepGadget()
  {
    if (host_data_buffer_) delete [] host_data_buffer_;
  }

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

    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    if( cfg.get() == 0x0 ){
      GADGET_DEBUG1("Unable to parse Ismrmrd header\n");
      return GADGET_FAIL;
    }

    use_multiframe_grouping_ = get_bool_value(std::string("use_multiframe_grouping").c_str());

    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();

    if (e_seq.size() != 1) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
      GADGET_DEBUG1("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    //ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

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

    image_dimensions_.push_back(r_space.matrixSize().x());
    image_dimensions_.push_back(r_space.matrixSize().y());

    fov_vec_.push_back(r_space.fieldOfView_mm().x());
    fov_vec_.push_back(r_space.fieldOfView_mm().y());
    fov_vec_.push_back(r_space.fieldOfView_mm().z());

    slices_ = e_limits.slice().present() ? e_limits.slice().get().maximum() + 1 : 1;
    sets_ = e_limits.set().present() ? e_limits.set().get().maximum() + 1 : 1;

    buffer_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);

    GADGET_DEBUG2("smax:                    %f\n", smax_);
    GADGET_DEBUG2("gmax:                    %f\n", gmax_);
    GADGET_DEBUG2("Tsamp_ns:                %d\n", Tsamp_ns_);
    GADGET_DEBUG2("Nints:                   %d\n", Nints_);
    GADGET_DEBUG2("fov:                     %f\n", fov_);
    GADGET_DEBUG2("krmax:                   %f\n", krmax_);
    GADGET_DEBUG2("samples_to_skip_start_ : %d\n", samples_to_skip_start_);
    GADGET_DEBUG2("samples_to_skip_end_   : %d\n", samples_to_skip_end_);
    GADGET_DEBUG2("matrix_size_x          : %d\n", image_dimensions_[0]);
    GADGET_DEBUG2("matrix_size_y          : %d\n", image_dimensions_[1]);

    return GADGET_OK;
  }

  int gpuSpiralSensePrepGadget::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    if (is_noise) { //Noise should have been consumed by the noise adjust, but just in case.
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

      try{host_traj_->create(&trajectory_dimensions);}
      catch (std::runtime_error &err){
	GADGET_DEBUG_EXCEPTION(err,"Unable to allocate memory for trajectory\n");
	return GADGET_FAIL;
      }

      try{host_weights_->create(&trajectory_dimensions);}
      catch (std::runtime_error& err ){
	GADGET_DEBUG_EXCEPTION(err,"Unable to allocate memory for weights\n");
	return GADGET_FAIL;
      }

      float* co_ptr = reinterpret_cast<float*>(host_traj_->get_data_ptr());
      float* we_ptr =  reinterpret_cast<float*>(host_weights_->get_data_ptr());

      for (int i = 0; i < (samples_per_interleave_*Nints_); i++) {
	co_ptr[i*2]   = -x_trajectory[i]/2;
	co_ptr[i*2+1] = -y_trajectory[i]/2;
	we_ptr[i] = weighting[i];
      }

      delete [] xgrad;
      delete [] ygrad;
      delete [] x_trajectory;
      delete [] y_trajectory;
      delete [] weighting;

      //Make NFFT plan
      // Matrix sizes
      uintd2 matrix_size = uintd2(image_dimensions_[0],image_dimensions_[1]);
      uintd2 matrix_size_os = uintd2(image_dimensions_[0]*2,image_dimensions_[1]*2);

      // Kernel width
      float W = 5.5f;

      // Upload host arrays to device arrays
      cuNDArray<floatd2> traj;
      try {traj= cuNDArray<floatd2>(*host_traj_);}
      catch (std::runtime_error& err){
	GADGET_DEBUG_EXCEPTION(err,"Failed to allocate device array\n");
	return GADGET_FAIL;
      }

      try{gpu_weights_ = cuNDArray<float>(*host_weights_);}
      catch (std::runtime_error& err){
	GADGET_DEBUG_EXCEPTION(err,"Failed to allocate device array\n");
	return GADGET_FAIL;
      };

      // Initialize plan
      // NFFT_plan<float, 2> plan( matrix_size, matrix_size_os, W );
      plan_ = cuNFFT_plan<float, 2>( matrix_size, matrix_size_os, W );

      // Preprocess
      try{plan_.preprocess( &traj, cuNFFT_plan<float,2>::NFFT_PREP_ALL );}
      catch (std::runtime_error& err){
	GADGET_DEBUG_EXCEPTION(err,"NFFT preprocess failed\n");
	return GADGET_FAIL;
      }

      prepared_ = true;
    }

    if (!host_data_buffer_) {
      std::vector<unsigned int> data_dimensions;
      data_dimensions.push_back(samples_per_interleave_*interleaves_);
      data_dimensions.push_back(m1->getObjectPtr()->active_channels);

      host_data_buffer_ = new hoNDArray<float_complext>[slices_*sets_];
      if (!host_data_buffer_) {
	GADGET_DEBUG1("Unable to allocate array for host data buffer\n");
	return GADGET_FAIL;
      }

      for (unsigned int i = 0; i < slices_*sets_; i++) {
	try{ host_data_buffer_[i].create(&data_dimensions); }
	catch (std::exception & err) {
	  GADGET_DEBUG1("Unable to allocate memory for data buffer\n");
	  return GADGET_FAIL;
	}
	host_data_buffer_[i].fill(0.0f);
      }
    }

    interleaves_counter_singleframe_++;
    interleaves_counter_multiframe_++;
    unsigned int samples_to_copy = m1->getObjectPtr()->number_of_samples-samples_to_skip_end_;
    unsigned int interleave = m1->getObjectPtr()->idx.kspace_encode_step_1;
    unsigned int slice = m1->getObjectPtr()->idx.slice;
    unsigned int set = m1->getObjectPtr()->idx.set;

    unsigned int samples_per_channel =  host_data_buffer_->get_size(0);

    buffer_[set*slices_+slice].enqueue_tail(m1);
    ISMRMRD::AcquisitionHeader base_head = *m1->getObjectPtr();

    if (samples_to_skip_end_ == -1) {
      samples_to_skip_end_ = m1->getObjectPtr()->number_of_samples-samples_per_interleave_;
      GADGET_DEBUG2("Adjusting samples_to_skip_end_ = %d\n", samples_to_skip_end_);
    }

    std::complex<float>* data_ptr    = reinterpret_cast< std::complex<float>* >(host_data_buffer_[set*slices_+slice].get_data_ptr());
    std::complex<float>* profile_ptr = m2->getObjectPtr()->get_data_ptr();

    for (unsigned int c = 0; c < m1->getObjectPtr()->active_channels; c++) {
      memcpy(data_ptr+c*samples_per_channel+interleave*samples_to_copy,
	     profile_ptr+c*m1->getObjectPtr()->number_of_samples, samples_to_copy*sizeof(std::complex<float>));
    }

    bool is_last_scan_in_slice = ISMRMRD::FlagBit(ISMRMRD::ACQ_LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags);
    if (is_last_scan_in_slice) {

      if( Nints_%interleaves_counter_singleframe_ ){
	GADGET_DEBUG1("Unexpected number of interleaves encountered in frame\n");
	return GADGET_FAIL;
      }

      if( acceleration_factor_ != Nints_/interleaves_counter_singleframe_ ){
	GADGET_DEBUG1("Change of acceleration factor detected\n");
	acceleration_factor_ =  Nints_/interleaves_counter_singleframe_;
      }

      if( !use_multiframe_grouping_ || (use_multiframe_grouping_ && interleaves_counter_multiframe_ == Nints_) ){

	GPUTimer timer("Spiral gridding for csm calc...");

	unsigned int num_coils = m1->getObjectPtr()->active_channels;

	cuNDArray<float_complext> data(&host_data_buffer_[set*slices_+slice]);

	// Setup averaged image (an image for each coil)
	std::vector<unsigned int> image_dims;
	image_dims.push_back(image_dimensions_[0]);
	image_dims.push_back(image_dimensions_[1]);
	image_dims.push_back(num_coils);
	cuNDArray<float_complext> image; 
	
	try{image.create(&image_dims);}
	catch (std::runtime_error &err){
	  GADGET_DEBUG_EXCEPTION(err,"\nError allocating coil images on device\n");
	  return GADGET_FAIL;
	}

	try {plan_.compute( &data, &image, &gpu_weights_, cuNFFT_plan<float,2>::NFFT_BACKWARDS_NC2C );}
	catch (std::runtime_error& err){
	  GADGET_DEBUG_EXCEPTION(err,"NFFT compute failed\n");
	  return GADGET_FAIL;
	}

	// Estimate CSM
	boost::shared_ptr< cuNDArray<float_complext> > csm = estimate_b1_map<float,2>( &image );

	// Use a SENSE operator to calculate a combined image using the trajectories.
	boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E(new cuNonCartesianSenseOperator<float,2>());
	E->set_csm(csm);

	// Setup averaged image (one image from the combined coils)
	boost::shared_ptr< std::vector<unsigned int> > reg_dims = image.get_dimensions();
	reg_dims->pop_back();
	cuNDArray<float_complext> reg_image;

	try{reg_image.create(reg_dims.get());}
	catch (std::runtime_error &err){
	  GADGET_DEBUG_EXCEPTION(err,"\nError allocating regularization image on device\n");
	  return GADGET_FAIL;
	}

	try {E->mult_csm_conj_sum( &image, &reg_image ); }
	catch (std::runtime_error& err){
	  GADGET_DEBUG_EXCEPTION(err,"\nError combining coils to regularization image\n");
	  return GADGET_FAIL;
	}

	boost::shared_ptr< hoNDArray<float_complext> > csm_host = csm->to_host();
	boost::shared_ptr< hoNDArray<float_complext> > reg_host = reg_image.to_host();

	unsigned int profiles_buffered = buffer_[set*slices_+slice].message_count();

	std::vector<unsigned int> ddimensions;

	ddimensions.push_back(samples_per_interleave_*interleaves_counter_singleframe_*((use_multiframe_grouping_) ? acceleration_factor_ : 1));
	ddimensions.push_back(num_coils);
	
	boost::shared_ptr< hoNDArray<float_complext> > data_host(new hoNDArray<float_complext>());

	try{data_host->create(&ddimensions);}
	catch (std::runtime_error& err){
	  GADGET_DEBUG_EXCEPTION(err,"Unable to allocate host data array\n");
	  return GADGET_FAIL;
	}

	ddimensions.clear();
	ddimensions.push_back(samples_per_interleave_*interleaves_counter_singleframe_);
	ddimensions.push_back((use_multiframe_grouping_) ? acceleration_factor_ : 1);

	boost::shared_ptr< hoNDArray<floatd2> > traj_host(new hoNDArray<floatd2>());
	try {traj_host->create(&ddimensions);}
	catch (std::runtime_error& err){
	  GADGET_DEBUG_EXCEPTION(err, "Unable to allocate host trajectory array\n");
	  return GADGET_FAIL;
	}

	boost::shared_ptr< hoNDArray<float> > dcw_host(new hoNDArray<float>());
	try {dcw_host->create(&ddimensions);}
	catch (std::runtime_error& err){
	  GADGET_DEBUG_EXCEPTION(err, "Unable to allocate host density compensation array\n");
	  return GADGET_FAIL;
	}

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

	if (buffer_[set*slices_+slice].message_count()) {
	  GADGET_DEBUG1("Error occured, all messages should have been cleared off the buffer by now.\n");
	  return GADGET_FAIL;
	}

	GadgetContainerMessage< SenseJob >* m4 = new GadgetContainerMessage< SenseJob >();

	m4->getObjectPtr()->dat_host_ = data_host;
	m4->getObjectPtr()->csm_host_ = csm_host;
	m4->getObjectPtr()->reg_host_ = reg_host;
	m4->getObjectPtr()->tra_host_ = traj_host;
	m4->getObjectPtr()->dcw_host_ = dcw_host;

	GadgetContainerMessage<ISMRMRD::ImageHeader>* m3 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	m3->cont(m4);

	m3->getObjectPtr()->matrix_size[0] = image_dimensions_[0];
	m3->getObjectPtr()->matrix_size[1] = image_dimensions_[1];
	m3->getObjectPtr()->matrix_size[2] = acceleration_factor_;

	m3->getObjectPtr()->field_of_view[0] = fov_vec_[0];
	m3->getObjectPtr()->field_of_view[1] = fov_vec_[1];
	m3->getObjectPtr()->field_of_view[2] = fov_vec_[2];

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

	if (this->next()->putq(m3) < 0) {
	  GADGET_DEBUG1("Failed to put job on queue.\n");
	  m3->release();
	  return GADGET_FAIL;
	}
	interleaves_counter_multiframe_ = 0;
      }
      interleaves_counter_singleframe_ = 0;
    }
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(gpuSpiralSensePrepGadget)
}
