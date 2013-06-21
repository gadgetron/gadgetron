#include "gpuCgSenseGadget.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "cuNDArray_utils.h"
#include "Gadgetron.h"
#include "GadgetMRIHeaders.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "vector_td_utilities.h"

//#include "hoNDArray_fileio.h"

namespace Gadgetron{

  gpuCgSenseGadget::gpuCgSenseGadget()
    : is_configured_(false)
    , channels_(0)
    , image_counter_(0)
  {
    set_parameter(std::string("deviceno").c_str(), "0");
    set_parameter(std::string("sliceno").c_str(), "0");
    set_parameter(std::string("number_of_iterations").c_str(), "5");
    set_parameter(std::string("cg_limit").c_str(), "1e-6");
    set_parameter(std::string("oversampling_factor").c_str(), "1.25");
    set_parameter(std::string("kernel_width").c_str(), "5.5");
    set_parameter(std::string("kappa").c_str(), "0.3");
    
    matrix_size_ = uintd2(0,0);
    matrix_size_os_ = uintd2(0,0);
    matrix_size_seq_ = uintd2(0,0);
  }

  gpuCgSenseGadget::~gpuCgSenseGadget() {}

  int gpuCgSenseGadget::process_config( ACE_Message_Block* mb )
  {
    //GADGET_DEBUG1("gpuCgSenseGadget::process_config\n");

    device_number_ = get_int_value(std::string("deviceno").c_str());

    int number_of_devices = 0;
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

    pass_on_undesired_data_ = get_bool_value(std::string("pass_on_undesired_data").c_str());
    slice_number_ = get_int_value(std::string("sliceno").c_str());
    number_of_iterations_ = get_int_value(std::string("number_of_iterations").c_str());
    cg_limit_ = get_double_value(std::string("cg_limit").c_str());
    oversampling_factor_ = get_double_value(std::string("oversampling_factor").c_str());
    kernel_width_ = get_double_value(std::string("kernel_width").c_str());
    kappa_ = get_double_value(std::string("kappa").c_str());

    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    std::vector<long> dims;
    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
    if (e_seq.size() != 1) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
      GADGET_DEBUG1("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    //ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    //ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    matrix_size_seq_ = uintd2( r_space.matrixSize().x(), r_space.matrixSize().y() );

    if (!is_configured_) {

      channels_ = cfg->acquisitionSystemInformation().present() ?
	(cfg->acquisitionSystemInformation().get().receiverChannels().present() ? cfg->acquisitionSystemInformation().get().receiverChannels().get() : 1) : 1;

      // Allocate encoding operator for non-Cartesian Sense
      E_ = boost::shared_ptr< cuNonCartesianSenseOperator<float,2> >( new cuNonCartesianSenseOperator<float,2>() );

      // Allocate preconditioner
      D_ = boost::shared_ptr< cuCgPreconditioner<float_complext> >( new cuCgPreconditioner<float_complext>() );

      // Allocate regularization image operator
      R_ = boost::shared_ptr< cuImageOperator<float_complext> >( new cuImageOperator<float_complext>() );
      R_->set_weight( kappa_ );

      // Setup solver
      cg_.set_encoding_operator( E_ );        // encoding matrix
      cg_.add_regularization_operator( R_ );  // regularization matrix
      cg_.set_preconditioner( D_ );           // preconditioning matrix
      cg_.set_max_iterations( number_of_iterations_ );
      cg_.set_tc_tolerance( cg_limit_ );
      cg_.set_output_mode( cuCgSolver<float_complext>::OUTPUT_SILENT);

      is_configured_ = true;
    }

    return GADGET_OK;
  }

  int gpuCgSenseGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<SenseJob> *m2)
  {
    // Is this data for this gadget's slice?
    //
    
    if (m1->getObjectPtr()->slice != slice_number_) {      
      // No, pass it downstream...
      return this->next()->putq(m1);
    }
    
    //GADGET_DEBUG1("gpuCgSenseGadget::process\n");
    //GPUTimer timer("gpuCgSenseGadget::process");

    if (!is_configured_) {
      GADGET_DEBUG1("Data received before configuration was completed\n");
      return GADGET_FAIL;
    }

    SenseJob* j = m2->getObjectPtr();

    // Some basic validation of the incoming Sense job
    if (!j->csm_host_.get() || !j->dat_host_.get() || !j->tra_host_.get() || !j->dcw_host_.get()) {
      GADGET_DEBUG1("Received an incomplete Sense job\n");
      return GADGET_FAIL;
    }

    unsigned int samples = j->dat_host_->get_size(0);
    unsigned int channels = j->dat_host_->get_size(1);
    unsigned int rotations = samples / j->tra_host_->get_number_of_elements();
    unsigned int frames = j->tra_host_->get_size(1)*rotations;

    if( samples%j->tra_host_->get_number_of_elements() ) {
      GADGET_DEBUG2("Mismatch between number of samples (%d) and number of k-space coordinates (%d).\nThe first should be a multiplum of the latter.\n", 
		    samples, j->tra_host_->get_number_of_elements());
      return GADGET_FAIL;
    }

    boost::shared_ptr< cuNDArray<floatd2> > traj(new cuNDArray<floatd2> (j->tra_host_.get()));
    boost::shared_ptr< cuNDArray<float> > dcw(new cuNDArray<float> (j->dcw_host_.get()));
    boost::shared_ptr< cuNDArray<float_complext> > csm(new cuNDArray<float_complext> (j->csm_host_.get()));
    boost::shared_ptr< cuNDArray<float_complext> > device_samples(new cuNDArray<float_complext> (j->dat_host_.get()));

    cudaDeviceProp deviceProp;
    if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
      GADGET_DEBUG1( "Error: unable to query device properties.\n" );
      return GADGET_FAIL;
    }
    
    unsigned int warp_size = deviceProp.warpSize;
    
    matrix_size_ = uintd2( j->reg_host_->get_size(0), j->reg_host_->get_size(1) );    

    matrix_size_os_ =
      uintd2(((static_cast<unsigned int>(std::ceil(matrix_size_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
	     ((static_cast<unsigned int>(std::ceil(matrix_size_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);
    
    GADGET_DEBUG2("Matrix size    : [%d,%d] \n", matrix_size_[0], matrix_size_[1]);    
    GADGET_DEBUG2("Matrix size OS : [%d,%d] \n", matrix_size_os_[0], matrix_size_os_[1]);

    std::vector<unsigned int> image_dims = to_std_vector(matrix_size_);
    image_dims.push_back(frames);
    
    E_->set_domain_dimensions(&image_dims);
    E_->set_codomain_dimensions(device_samples->get_dimensions().get());
    E_->set_dcw(dcw);
    E_->set_csm(csm);

    E_->setup( matrix_size_, matrix_size_os_, static_cast<float>(kernel_width_) );
    E_->preprocess(traj.get());

    boost::shared_ptr< cuNDArray<float_complext> > reg_image(new cuNDArray<float_complext> (j->reg_host_.get()));
    R_->compute(reg_image.get());

    // Define preconditioning weights
    boost::shared_ptr< cuNDArray<float> > _precon_weights = sum(abs_square(csm.get()).get(), 2);
    boost::shared_ptr<cuNDArray<float> > R_diag = R_->get();
    *R_diag *= float(kappa_);
    *_precon_weights += *R_diag;
    R_diag.reset();
    reciprocal_sqrt_inplace(_precon_weights.get());	
    boost::shared_ptr< cuNDArray<float_complext> > precon_weights = real_to_complex<float_complext>( _precon_weights.get() );
    _precon_weights.reset();
    D_->set_weights( precon_weights );
	
    // Invoke solver
    // 

    boost::shared_ptr< cuNDArray<float_complext> > cgresult;
    
    {
      //GPUTimer timer("gpuCgSenseGadget::solve()");
      cgresult = cg_.solve(device_samples.get());
    }

    if (!cgresult.get()) {
      GADGET_DEBUG1("Iterative_sense_compute failed\n");
      return GADGET_FAIL;
    }

    /*
    static int counter = 0;
    char filename[256];
    sprintf((char*)filename, "recon_%d.real", counter);
    write_nd_array<float>( abs(cgresult.get())->to_host().get(), filename );
    counter++; */

    // If the recon matrix size exceeds the sequence matrix size then crop
    if( matrix_size_seq_ != matrix_size_ )
      cgresult = crop<float_complext,2>( (matrix_size_-matrix_size_seq_)>>2, matrix_size_seq_, cgresult.get() );    
    
    // Now pass on the reconstructed images

    for( unsigned int frame=0; frame<frames; frame++ ){
      
      GadgetContainerMessage<ISMRMRD::ImageHeader> *m = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
      GadgetContainerMessage< hoNDArray< std::complex<float> > > *cm = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();      
      
      if( !m || !cm ){
	GADGET_DEBUG1("Unable create container messages\n");
	return GADGET_FAIL;
      }

      *(m->getObjectPtr()) = *(m1->getObjectPtr());
      m->cont(cm);
      
      std::vector<unsigned int> img_dims(2);
      img_dims[0] = matrix_size_seq_[0];
      img_dims[1] = matrix_size_seq_[1];

      cm->getObjectPtr()->create(&img_dims);

      size_t data_length = prod(matrix_size_seq_);

      cudaMemcpy(cm->getObjectPtr()->get_data_ptr(),
		 cgresult->get_data_ptr()+frame*data_length,
		 data_length*sizeof(std::complex<float>),
		 cudaMemcpyDeviceToHost);

      cudaError_t err = cudaGetLastError();
      if( err != cudaSuccess ){
	GADGET_DEBUG2("Unable to copy result from device to host: %s\n", cudaGetErrorString(err));
	m->release();
	return GADGET_FAIL;
      }

      m->getObjectPtr()->matrix_size[0] = matrix_size_seq_[0];
      m->getObjectPtr()->matrix_size[1] = matrix_size_seq_[1];
      m->getObjectPtr()->matrix_size[2] = 1;
      m->getObjectPtr()->channels       = 1;
      
      if (this->next()->putq(m) < 0) {
	GADGET_DEBUG1("Failed to put result image on to queue\n");
	m->release();
	return GADGET_FAIL;
      }
    }
    
    m1->release();
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(gpuCgSenseGadget)
}
