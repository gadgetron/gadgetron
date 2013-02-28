#include "cuNDArray.h"
#include "Gadgetron.h"
#include "SpiralGadget.h"
#include "GadgetXml.h"
#include "hoNDArray_fileio.h"
#include "vector_td.h"
#include "vector_td_operators.h"
#include "check_CUDA.h"

#include <vector>
namespace Gadgetron{
void calc_vds(double slewmax,double gradmax,double Tgsample,double Tdsample,int Ninterleaves,
		double* fov, int numfov,double krmax,
		int ngmax, double** xgrad,double** ygrad,int* numgrad);

void calc_traj(double* xgrad, double* ygrad, int ngrad, int Nints, double Tgsamp, double krmax,
		double** x_trajectory, double** y_trajectory,
		double** weights);


SpiralGadget::SpiralGadget()
: samples_to_skip_start_(0)
, samples_to_skip_end_(0)
, samples_per_interleave_(0)
, host_data_buffer_(0)
, image_counter_(0)
, image_series_(0)
{
	GADGET_DEBUG1("Initializing Spiral\n");
}

SpiralGadget::~SpiralGadget()
{
	if (host_data_buffer_) delete [] host_data_buffer_;
}

int SpiralGadget::process_config(ACE_Message_Block* mb)
{
	TiXmlDocument doc;
	doc.Parse(mb->rd_ptr());

	GADGET_DEBUG1("Calculating trajectory\n");

	GadgetXMLNode n = GadgetXMLNode(&doc).get<GadgetXMLNode>(std::string("gadgetron"))[0];

	int     Tsamp_ns = n.get<long>(std::string("wip.long.value"))[4];		//"samplingtime_ns.value"))[0];
	int     Nints  = n.get<long>(std::string("wip.long.value"))[11];		//"interleaves.value"))[0];
	double  gmax   = n.get<double>(std::string("wip.double.value"))[0];		//"maxgradient_gcms.value"))[0];
	double  smax   = n.get<double>(std::string("wip.double.value"))[1];		//"maxslewrate_gcmsi.value"))[0];
	double  krmax  = n.get<double>(std::string("wip.double.value"))[3];		//"krmax_cm.value"))[0];
	double  fov    = n.get<double>(std::string("wip.double.value"))[4];		//"fovcoeff_1.value"))[0];
	int     nfov   = 1;         /*  number of fov coefficients.             */
	int     ngmax  = 1e5;       /*  maximum number of gradient samples      */
	double  *xgrad;             /*  x-component of gradient.                */
	double  *ygrad;             /*  y-component of gradient.                */
	double  *x_trajectory;
	double  *y_trajectory;
	double  *weighting;
	int     ngrad;
	//int     count;
	double sample_time = (1.0*Tsamp_ns) * 1e-9;


	samples_to_skip_start_  = 0; //n.get<int>(std::string("samplestoskipstart.value"))[0];
	samples_to_skip_end_    = -1; //n.get<int>(std::string("samplestoskipend.value"))[0];
	//  samples_per_adc_        = 0; //n.get<int>(std::string("samplesperadc.value"))[0];
	//  adcs_per_interleave_    = 1; //n.get<int>(std::string("adcsperinterleave.value"))[0];

	image_dimensions_.push_back(n.get<long>(std::string("encoding.kspace.matrix_size.value"))[0]);
	image_dimensions_.push_back(n.get<long>(std::string("encoding.kspace.matrix_size.value"))[1]);

	GADGET_DEBUG2("smax:                    %f\n", smax);
	GADGET_DEBUG2("gmax:                    %f\n", gmax);
	GADGET_DEBUG2("Tsamp_ns:                %d\n", Tsamp_ns);
	GADGET_DEBUG2("sample_time:             %f\n", sample_time);
	GADGET_DEBUG2("Nints:                   %d\n", Nints);
	GADGET_DEBUG2("fov:                     %f\n", fov);
	GADGET_DEBUG2("krmax:                   %f\n", krmax);
	GADGET_DEBUG2("samples_to_skip_start_ : %d\n", samples_to_skip_start_);
	GADGET_DEBUG2("samples_to_skip_end_   : %d\n", samples_to_skip_end_);
	GADGET_DEBUG2("matrix_size_x          : %d\n", image_dimensions_[0]);
	GADGET_DEBUG2("matrix_size_y          : %d\n", image_dimensions_[1]);



	/*	call c-function here to calculate gradients */
	calc_vds(smax,gmax,sample_time,sample_time,Nints,&fov,nfov,krmax,ngmax,&xgrad,&ygrad,&ngrad);
	samples_per_interleave_ = ngrad;


	/* Calcualte the trajectory and weights*/
	calc_traj(xgrad, ygrad, ngrad, Nints, sample_time, krmax, &x_trajectory, &y_trajectory, &weighting);

	host_traj_ = boost::shared_ptr< hoNDArray<floatd2> >(new hoNDArray<floatd2>);
	host_weights_ = boost::shared_ptr< hoNDArray<float> >(new hoNDArray<float>);

	std::vector<unsigned int> trajectory_dimensions;
	trajectory_dimensions.push_back(samples_per_interleave_*Nints);

	try{host_traj_->create(&trajectory_dimensions);}
	catch (runtime_error &err ){
		GADGET_DEBUG_EXCEPTION(err,"Unable to allocate memory for trajectory\n");
		return GADGET_FAIL;
	}

	try{host_weights_->create(&trajectory_dimensions);}
	catch (runtime_error &err ){
		GADGET_DEBUG_EXCEPTION(err,"Unable to allocate memory for weights\n");
		return GADGET_FAIL;
	}


	float* co_ptr = reinterpret_cast<float*>(host_traj_->get_data_ptr());
	float* we_ptr =  reinterpret_cast<float*>(host_weights_->get_data_ptr());

	for (int i = 0; i < (ngrad*Nints); i++) {
		co_ptr[i*2]   = -x_trajectory[i]/2;
		co_ptr[i*2+1] = -y_trajectory[i]/2;
		we_ptr[i] = weighting[i];
	}

	delete [] xgrad;
	delete [] ygrad;
	delete [] x_trajectory;
	delete [] y_trajectory;
	delete [] weighting;

	unsigned int slices = n.get<long>(std::string("encoding.slices.value"))[0];

	std::vector<unsigned int> data_dimensions;
	data_dimensions.push_back(ngrad*Nints);
	data_dimensions.push_back(n.get<long>(std::string("encoding.channels.value"))[0]);

	host_data_buffer_ = new hoNDArray<float_complext>[slices];
	if (!host_data_buffer_) {
		GADGET_DEBUG1("Unable to allocate array for host data buffer\n");
		return GADGET_FAIL;
	}

	for (unsigned int i = 0; i < slices; i++) {
		try{host_data_buffer_[i].create(&data_dimensions);}
		catch (runtime_error &err){
			GADGET_DEBUG_EXCEPTION(err,"Unable to allocate memory for data buffer\n");
			return GADGET_FAIL;
		}
	}


	//Make NFFT plan
	// Matrix sizes
	uintd2 matrix_size = uintd2(image_dimensions_[0],image_dimensions_[1]);
	uintd2 matrix_size_os = uintd2(image_dimensions_[0]*2,image_dimensions_[1]*2);

	// Kernel width
	float W = 5.5f;

	// Upload host arrays to device arrays
	cuNDArray<floatd2> traj(host_traj_.get());
	gpu_weights_ = cuNDArray<float>(host_weights_.get());

	// Initialize plan
	// NFFT_plan<float, 2> plan( matrix_size, matrix_size_os, W );
	plan_ = NFFT_plan<float, 2>( matrix_size, matrix_size_os, W );

	// Preprocess
	try { 	plan_.preprocess( &traj, NFFT_plan<float,2>::NFFT_PREP_ALL ); }
	catch (runtime_error& err){
		GADGET_DEBUG_EXCEPTION(err,"NFFT preprocess failed\n");
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

int SpiralGadget::
process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

	if (samples_to_skip_end_ == -1) {
		samples_to_skip_end_ = m1->getObjectPtr()->number_of_samples-samples_per_interleave_;
		GADGET_DEBUG2("Adjusting samples_to_skip_end_ = %d\n", samples_to_skip_end_);
	}

	unsigned int samples_to_copy = m1->getObjectPtr()->number_of_samples-samples_to_skip_end_;
	unsigned int interleave = m1->getObjectPtr()->idx.kspace_encode_step_1;
	unsigned int slice = m1->getObjectPtr()->idx.slice;

	unsigned int samples_per_channel =  host_data_buffer_->get_size(0);

	std::complex<float>* data_ptr    = reinterpret_cast< std::complex<float>* >(host_data_buffer_[slice].get_data_ptr());
	std::complex<float>* profile_ptr = m2->getObjectPtr()->get_data_ptr();

	for (unsigned int c = 0; c < m1->getObjectPtr()->active_channels; c++) {
		memcpy(data_ptr+c*samples_per_channel+interleave*samples_to_copy,
				profile_ptr+c*m1->getObjectPtr()->number_of_samples, samples_to_copy*sizeof(std::complex<float>));
	}

	if (ISMRMRD::FlagBit(ISMRMRD::ACQ_LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags)) {

		unsigned int num_batches = m1->getObjectPtr()->active_channels;

		cuNDArray<float_complext> data(&host_data_buffer_[slice]);

		// Setup result image
		std::vector<unsigned int> image_dims;
		image_dims.push_back(image_dimensions_[0]);
		image_dims.push_back(image_dimensions_[1]);
		image_dims.push_back(num_batches);
		cuNDArray<float_complext> image; image.create(&image_dims);

		try{ plan_.compute( &data, &image, &gpu_weights_, NFFT_plan<float,2>::NFFT_BACKWARDS_NC2C ); }
		catch (runtime_error& err){
			GADGET_DEBUG_EXCEPTION(err, "NFFT compute failed\n");
			return GADGET_FAIL;
		}

		boost::shared_ptr< hoNDArray<float_complext> > image_host = image.to_host();

		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m4 =
				new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

		try{ m4->getObjectPtr()->create(&image_dimensions_);}
		catch (runtime_error& err){
			GADGET_DEBUG_EXCEPTION(err,"Unable to allocate memory for combined image\n");
			m4->release();
			return GADGET_FAIL;
		}

		unsigned int npixels = image_dimensions_[0]*image_dimensions_[1];
		std::complex<float>* recon_ptr    = reinterpret_cast< std::complex<float>* >(image_host->get_data_ptr());
		std::complex<float>* comb_ptr     = reinterpret_cast< std::complex<float>* >(m4->getObjectPtr()->get_data_ptr());

		for (unsigned int i = 0; i < npixels; i++) {
			float mag = 0.0;
			float phase = 0.0;
			for (unsigned int c = 0; c < num_batches; c++) {
				float mag_tmp = norm(recon_ptr[c*npixels+i]);
				phase += mag_tmp*arg(recon_ptr[c*npixels+i]);
				mag += mag_tmp;
			}
			comb_ptr[i] = std::polar(std::sqrt(mag),phase)*std::complex<float>(npixels,0.0);
		}


		GadgetContainerMessage<ISMRMRD::ImageHeader>* m3 =
				new GadgetContainerMessage<ISMRMRD::ImageHeader>();

		m3->cont(m4);

		m3->getObjectPtr()->matrix_size[0] = image_dimensions_[0];
		m3->getObjectPtr()->matrix_size[1] = image_dimensions_[1];
		m3->getObjectPtr()->matrix_size[2] = 1;
		m3->getObjectPtr()->channels       = 1;
		m3->getObjectPtr()->slice          = m1->getObjectPtr()->idx.slice;

		memcpy(m3->getObjectPtr()->position,m1->getObjectPtr()->position,
				sizeof(float)*3);

		memcpy(m3->getObjectPtr()->quaternion,m1->getObjectPtr()->quaternion,
				sizeof(float)*4);

		memcpy(m3->getObjectPtr()->patient_table_position, m1->getObjectPtr()->patient_table_position, sizeof(float)*3);

		m3->getObjectPtr()->image_data_type = ISMRMRD::DATA_COMPLEX_FLOAT;
		m3->getObjectPtr()->image_index = ++image_counter_; 
		m3->getObjectPtr()->image_series_index = image_series_;

		if (this->next()->putq(m3) < 0) {
			m3->release();
			return GADGET_FAIL;
		}

	}
	m1->release();
	return GADGET_OK;
}


GADGET_FACTORY_DECLARE(SpiralGadget)
}
