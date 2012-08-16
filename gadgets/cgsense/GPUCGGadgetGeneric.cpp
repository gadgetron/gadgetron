#include "GPUCGGadgetGeneric.h"
#include "Gadgetron.h"
#include "GadgetMRIHeaders.h"
#include "ndarray_vector_td_utilities.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "GadgetIsmrmrdReadWrite.h"

#include "hoNDArray_fileio.h"

#include "tinyxml.h"

GPUCGGadgetGeneric::GPUCGGadgetGeneric()
: channels_(0)
, device_number_(0)
, number_of_iterations_(5)
, cg_limit_(1e-6)
, oversampling_(1.25)
, kernel_width_(5.5)
, kappa_(0.1)
, is_configured_(false)
, image_series_(0)
, image_counter_(0)
{
	matrix_size_ = uintd2(0,0);
	matrix_size_os_ = uintd2(0,0);
}

GPUCGGadgetGeneric::~GPUCGGadgetGeneric() {}

int GPUCGGadgetGeneric::process_config( ACE_Message_Block* mb )
{
	GADGET_DEBUG1("GPUCGGadgetGeneric::process_config\n");


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

	number_of_iterations_ = get_int_value(std::string("number_of_iterations").c_str());
	cg_limit_ = get_double_value(std::string("cg_limit").c_str());
	oversampling_ = get_double_value(std::string("oversampling").c_str());
	kernel_width_ = get_double_value(std::string("kernel_width").c_str());
	kappa_ = get_double_value(std::string("kappa").c_str());
	pass_on_undesired_data_ = get_bool_value(std::string("pass_on_undesired_data").c_str());
	image_series_ = this->get_int_value("image_series");

	boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

	std::vector<long> dims;
	ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
	if (e_seq.size() != 1) {
		GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
		GADGET_DEBUG1("This Gadget only supports one encoding space\n");
		return GADGET_FAIL;
	}

	ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
	ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
	ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

	if (!is_configured_) {

		cudaDeviceProp deviceProp;
		if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
			GADGET_DEBUG1( "\nError: unable to query device properties.\n" );
			return GADGET_FAIL;
		}

		unsigned int warp_size = deviceProp.warpSize;

		channels_ = cfg->acquisitionSystemInformation().present() ?
					(cfg->acquisitionSystemInformation().get().receiverChannels().present() ? cfg->acquisitionSystemInformation().get().receiverChannels().get() : 1) : 1;

		matrix_size_ = uintd2(e_space.matrixSize().x(), e_space.matrixSize().y());

		GADGET_DEBUG2("Matrix size  : [%d,%d] \n", matrix_size_.vec[0], matrix_size_.vec[1]);

		matrix_size_os_ =
				uintd2(static_cast<unsigned int>(ceil((matrix_size_.vec[0]*oversampling_)/warp_size)*warp_size),
						static_cast<unsigned int>(ceil((matrix_size_.vec[1]*oversampling_)/warp_size)*warp_size));

		GADGET_DEBUG2("Matrix size OS: [%d,%d] \n", matrix_size_os_.vec[0], matrix_size_os_.vec[1]);

		// Allocate encoding operator for non-Cartesian Sense
		std::vector<unsigned int> image_dims = uintd_to_vector<2>(matrix_size_);
		E_ = boost::shared_ptr< cuNonCartesianSenseOperator<float,2> >( new cuNonCartesianSenseOperator<float,2>() );
		E_->set_device(device_number_);
		E_->set_domain_dimensions(&image_dims);

		// Allocate preconditioner
		D_ = boost::shared_ptr< cuCgPrecondWeights<float_complext> >( new cuCgPrecondWeights<float_complext>() );
		//D_->set_device(device_number_);

		// Allocate regularization image operator
		R_ = boost::shared_ptr< cuImageOperator<float,float_complext> >( new cuImageOperator<float,float_complext>() );
		R_->set_device(device_number_);
		R_->set_weight( kappa_ );
		cg_.set_device(device_number_);

		// Setup solver
		cg_.set_encoding_operator( E_ );        // encoding matrix
		cg_.add_regularization_operator( R_ );  // regularization matrix
		cg_.set_preconditioner( D_ );           // preconditioning matrix
		cg_.set_max_iterations( number_of_iterations_ );
		cg_.set_tc_tolerance( cg_limit_ );
		cg_.set_output_mode( cuCgSolver<float, float_complext>::OUTPUT_WARNINGS);

		if( configure_channels() == GADGET_FAIL )
			return GADGET_FAIL;

		is_configured_ = true;
	}

	return GADGET_OK;
}

int GPUCGGadgetGeneric::configure_channels()
{
	// We do not have a csm yet, so initialize a dummy one to purely ones
	boost::shared_ptr< cuNDArray<float_complext> > csm = boost::shared_ptr< cuNDArray<float_complext> >( new cuNDArray<float_complext> );
	std::vector<unsigned int> csm_dims = uintd_to_vector<2>(matrix_size_); csm_dims.push_back( channels_ );

	if( csm->create( &csm_dims ) == 0x0 ) {
		GADGET_DEBUG1( "\nError: unable to create csm.\n" );
		return GADGET_FAIL;
	}

	if( !cuNDA_clear<float_complext>( csm.get(), float_complext(1) ) ){
		GADGET_DEBUG1( "\nError: unable to clear csm.\n" );
		return GADGET_FAIL;
	}

	// Setup matrix operator
	E_->set_csm(csm);

	if( E_->setup( matrix_size_, matrix_size_os_, static_cast<float>(kernel_width_) ) < 0 ){
		GADGET_DEBUG1( "\nError: unable to setup encoding operator.\n" );
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

int GPUCGGadgetGeneric::process(GadgetContainerMessage<GadgetMessageImage>* m1, GadgetContainerMessage< CGSenseJob > * m2)
{

	//GPUTimer timer("GPUCGGadgetGeneric::process");

	if (!is_configured_) {
		GADGET_DEBUG1("\nData received before configuration complete\n");
		return GADGET_FAIL;
	}

	CGSenseJob* j = m2->getObjectPtr();

	//Let's first check that this job has the required stuff...
	if (!j->csm_host_.get() || !j->dat_host_.get() || !j->tra_host_.get() || !j->dcw_host_.get()) {
		GADGET_DEBUG1("Received an incomplete CGSense JOB\n");
		m1->release();
		return GADGET_FAIL;
	}

	unsigned int samples = j->dat_host_->get_size(0);
	unsigned int channels = j->dat_host_->get_size(1);

	if (samples != j->tra_host_->get_number_of_elements()) {
		GADGET_DEBUG2("Mismatch between number of samples (%d) and number of k-space coordinates (%d)\n", samples, j->tra_host_->get_number_of_elements());
		m1->release();
		return GADGET_FAIL;
	}


	if( m1->getObjectPtr()->channels != channels_ ) {
		GADGET_DEBUG2("Adjusting #channels from %d to %d\n", channels_,  m1->getObjectPtr()->channels );
		channels_ = m1->getObjectPtr()->channels;
		if( configure_channels() == GADGET_FAIL ) // Update buffers dependant on #channels
			return GADGET_FAIL;
	}

	boost::shared_ptr< cuNDArray<floatd2> > traj(new cuNDArray<floatd2> (j->tra_host_.get()));
	boost::shared_ptr< cuNDArray<float> > dcw(new cuNDArray<float> (j->dcw_host_.get()));
	boost::shared_ptr< cuNDArray<float_complext> > csm(new cuNDArray<float_complext> (j->csm_host_.get()));
	boost::shared_ptr< cuNDArray<float_complext> > device_samples(new cuNDArray<float_complext> (j->dat_host_.get()));

	E_->set_dcw(dcw);
	E_->set_csm(csm);
	if( E_->preprocess(traj.get()) < 0 ) {
		GADGET_DEBUG1("\nError during cgOperatorNonCartesianSense::preprocess()\n");
		return GADGET_FAIL;
	}


	boost::shared_ptr< cuNDArray<float_complext> > reg_image(new cuNDArray<float_complext> (j->reg_host_.get()));

	R_->compute(reg_image.get());

	// TODO: error check these computations

	// Define preconditioning weights
	boost::shared_ptr< cuNDArray<float> > _precon_weights = cuNDA_ss<float,float_complext>( csm.get(), 2 );
	cuNDA_axpy<float>( kappa_, R_->get(), _precon_weights.get() );
	cuNDA_reciprocal_sqrt<float>( _precon_weights.get() );
	boost::shared_ptr< cuNDArray<float_complext> > precon_weights = cuNDA_real_to_complext<float>( _precon_weights.get() );
	_precon_weights.reset();
	D_->set_weights( precon_weights );

	// Invoke solver
	boost::shared_ptr< cuNDArray<float_complext> > cgresult = cg_.solve(device_samples.get());

	if (!cgresult.get()) {
		GADGET_DEBUG1("\nIterative_sense_compute failed\n");
		return GADGET_FAIL;
	}

	m2->release();

	//Now pass the reconstructed image on

	GadgetContainerMessage< hoNDArray< std::complex<float> > >* cm2 =
			new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

	m1->cont(cm2);

	std::vector<unsigned int> img_dims(2);
	img_dims[0] = matrix_size_.vec[0];
	img_dims[1] = matrix_size_.vec[1];

	if (cm2->getObjectPtr()->create(&img_dims) == 0x0) {
		GADGET_DEBUG1("\nUnable to allocate host image array");
		m1->release();
		return GADGET_FAIL;
	}

	size_t data_length = prod(matrix_size_);

	cudaMemcpy(cm2->getObjectPtr()->get_data_ptr(),
			cgresult->get_data_ptr(),
			data_length*sizeof(std::complex<float>),
			cudaMemcpyDeviceToHost);

	cudaError_t err = cudaGetLastError();
	if( err != cudaSuccess ){
		GADGET_DEBUG2("\nUnable to copy result from device to host: %s", cudaGetErrorString(err));
		m1->release();
		return GADGET_FAIL;
	}

	m1->getObjectPtr()->matrix_size[0] = img_dims[0];
	m1->getObjectPtr()->matrix_size[1] = img_dims[1];
	m1->getObjectPtr()->matrix_size[2] = 1;
	m1->getObjectPtr()->channels       = 1;


	if (this->next()->putq(m1) < 0) {
		GADGET_DEBUG1("\nFailed to result image on to Q\n");
		m1->release();
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GPUCGGadgetGeneric)
