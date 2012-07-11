#include "Gadgetron.h"
#include "SpiralGadgetSW.h"
#include "GadgetXml.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "vector_td.h"
#include "vector_td_operators.h"
#include "check_CUDA.h"
#include "b1_map.h"
#include "cuNonCartesianSenseOperator.h"
#include "GPUCGGadgetGeneric.h"

#include <vector>

void calc_vds(double slewmax,double gradmax,double Tgsample,double Tdsample,int Ninterleaves,
		double* fov, int numfov,double krmax,
		int ngmax, double** xgrad,double** ygrad,int* numgrad);

void calc_traj(double* xgrad, double* ygrad, int ngrad, int Nints, double Tgsamp, double krmax,
		double** x_trajectory, double** y_trajectory,
		double** weights);


template <class T> int get_gadget_xml_parameter(GadgetXMLNode& n, std::string p, T& result, unsigned int index = 0)
{
	std::vector<T> parms = n.get<T>(p);
	if (parms.size() <= index) {
		GADGET_DEBUG2("Tried to get parameter %s from node and failed!.\n", p.c_str());
		return GADGET_FAIL;
	}

	result = parms[index];

	return GADGET_OK;
}

SpiralGadgetSW::SpiralGadgetSW()
: samples_to_skip_start_(0)
, samples_to_skip_end_(0)
, samples_per_interleave_(0)
, host_data_buffer_(0)
, image_counter_(0)
, image_series_(0)
{
	GADGET_DEBUG1("Initializing Spiral\n");
}

SpiralGadgetSW::~SpiralGadgetSW()
{
	if (host_data_buffer_) delete [] host_data_buffer_;
}

int SpiralGadgetSW::process_config(ACE_Message_Block* mb)
{

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


	TiXmlDocument doc;
	doc.Parse(mb->rd_ptr());


	GADGET_DEBUG1("Calculating trajectory\n");

	GadgetXMLNode n = GadgetXMLNode(&doc).get<GadgetXMLNode>(std::string("gadgetron"))[0];

	long     Tsamp_ns;
	if (get_gadget_xml_parameter(n,std::string("encoding.acquisition_dwell_time_ns.value"),Tsamp_ns) != GADGET_OK) {
		return GADGET_FAIL;
	}

	long     Nints;
	if (get_gadget_xml_parameter(n,std::string("spirals.Interleaves.value"),Nints) != GADGET_OK) {
		return GADGET_FAIL;
	}

	interleaves_ = static_cast<int>(Nints);

	double  gmax;
	if (get_gadget_xml_parameter(n,std::string("spirals.MaxGradient_Gcms.value"),gmax) != GADGET_OK) {
		return GADGET_FAIL;
	}

	double  smax;
	if (get_gadget_xml_parameter(n,std::string("spirals.MaxSlewRate_Gcmsi.value"),smax) != GADGET_OK) {
		return GADGET_FAIL;
	}

	double  krmax;
	if (get_gadget_xml_parameter(n,std::string("spirals.krmax_cm.value"),krmax) != GADGET_OK) {
		return GADGET_FAIL;
	}

	double  fov;
	if (get_gadget_xml_parameter(n,std::string("spirals.FOVCoeff_1.value"),fov) != GADGET_OK) {
		return GADGET_FAIL;
	}

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

	if (!host_traj_->create(&trajectory_dimensions)) {
		GADGET_DEBUG1("Unable to allocate memory for trajectory\n");
		return GADGET_FAIL;
	}

	if (!host_weights_->create(&trajectory_dimensions)) {
		GADGET_DEBUG1("Unable to allocate memory for weights\n");
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

	slices_ = static_cast<int>(n.get<long>(std::string("encoding.slices.value"))[0]);


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
	bool success = plan_.preprocess( &traj, NFFT_plan<float,2>::NFFT_PREP_ALL );

	if (!success) {
		GADGET_DEBUG1("NFFT preprocess failed\n");
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

int SpiralGadgetSW::
process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{


	if (!host_data_buffer_) {
		std::vector<unsigned int> data_dimensions;
		data_dimensions.push_back(samples_per_interleave_*interleaves_);
		data_dimensions.push_back(m1->getObjectPtr()->channels);

		host_data_buffer_ = new hoNDArray<float_complext>[slices_];
		if (!host_data_buffer_) {
			GADGET_DEBUG1("Unable to allocate array for host data buffer\n");
			return GADGET_FAIL;
		}

		for (unsigned int i = 0; i < slices_; i++) {
			if (!host_data_buffer_[i].create(&data_dimensions)) {
				GADGET_DEBUG1("Unable to allocate memory for data buffer\n");
				return GADGET_FAIL;
			}
		}
	}


	buffer_.enqueue_tail(m1);

	if (samples_to_skip_end_ == -1) {
		samples_to_skip_end_ = m1->getObjectPtr()->samples-samples_per_interleave_;
		GADGET_DEBUG2("Adjusting samples_to_skip_end_ = %d\n", samples_to_skip_end_);
	}

	unsigned int samples_to_copy = m1->getObjectPtr()->samples-samples_to_skip_end_;
	unsigned int interleave = m1->getObjectPtr()->idx.line;
	unsigned int slice = m1->getObjectPtr()->idx.slice;

	unsigned int samples_per_channel =  host_data_buffer_->get_size(0);

	std::complex<float>* data_ptr    = reinterpret_cast< std::complex<float>* >(host_data_buffer_[slice].get_data_ptr());
	std::complex<float>* profile_ptr = m2->getObjectPtr()->get_data_ptr();

	for (unsigned int c = 0; c < m1->getObjectPtr()->channels; c++) {
		memcpy(data_ptr+c*samples_per_channel+interleave*samples_to_copy,
				profile_ptr+c*m1->getObjectPtr()->samples, samples_to_copy*sizeof(std::complex<float>));
	}

	if (m1->getObjectPtr()->flags & GADGET_FLAG_LAST_ACQ_IN_SLICE) {

		unsigned int num_batches = m1->getObjectPtr()->channels;

		cuNDArray<float_complext> data(&host_data_buffer_[slice]);

		// Setup result image
		std::vector<unsigned int> image_dims;
		image_dims.push_back(image_dimensions_[0]);
		image_dims.push_back(image_dimensions_[1]);
		image_dims.push_back(num_batches);
		cuNDArray<float_complext> image; image.create(&image_dims);

		bool  success = plan_.compute( &data, &image, &gpu_weights_, NFFT_plan<float,2>::NFFT_BACKWARDS_NC2C );
		if (!success) {
			GADGET_DEBUG1("NFFT compute failed\n");
			return GADGET_FAIL;
		}


		// Estimate CSM
		boost::shared_ptr< cuNDArray<float_complext> > csm = estimate_b1_map<float,2>( &image );

		//We will use a SENSE operator to calculate a combined image using the trajectories.
		boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E_(new cuNonCartesianSenseOperator<float,2>());
		E_->set_csm(csm);

		boost::shared_ptr< std::vector<unsigned int> > reg_dims = image.get_dimensions();
		reg_dims->pop_back();

		cuNDArray<float_complext> reg_image;
		if( reg_image.create(reg_dims.get()) == 0x0 ){
			GADGET_DEBUG1("\nError allocating regularization image on device\n");
			return GADGET_FAIL;
		}

		if( E_->mult_csm_conj_sum( &image, &reg_image ) < 0 ){
			GADGET_DEBUG1("\nError combining coils to regularization image\n");
			return GADGET_FAIL;
		}

		boost::shared_ptr< hoNDArray<float_complext> > csm_host = csm->to_host();
		boost::shared_ptr< hoNDArray<float_complext> > reg_host = reg_image.to_host();


		unsigned int profiles_buffered = buffer_.message_count();
		boost::shared_ptr< hoNDArray<float_complext> > data_host(new hoNDArray<float_complext>());

		std::vector<unsigned int> ddimensions(2,0);
		ddimensions[0] = samples_per_interleave_*profiles_buffered;
		ddimensions[1] = num_batches; //Channels
		if (!data_host->create(&ddimensions)) {
			GADGET_DEBUG1("Unable to allocate host data array\n");
			return GADGET_FAIL;
		}

		ddimensions.pop_back();


		boost::shared_ptr< hoNDArray<floatd2> > traj_host(new hoNDArray<floatd2>());
		if (!traj_host->create(&ddimensions)) {
			GADGET_DEBUG1("Unable to allocate host trajectory array\n");
			return GADGET_FAIL;
		}

		boost::shared_ptr< hoNDArray<float> > dcw_host(new hoNDArray<float>());
		if (!dcw_host->create(&ddimensions)) {
			GADGET_DEBUG1("Unable to allocate host density compensation array\n");
			return GADGET_FAIL;
		}

		for (unsigned int p = 0; p < profiles_buffered; p++) {
			ACE_Message_Block* mbq;
			if (buffer_.dequeue_head(mbq) < 0) {
				GADGET_DEBUG1("Message dequeue failed\n");
				return GADGET_FAIL;
			}

			GadgetContainerMessage<GadgetMessageAcquisition>* acq = AsContainerMessage<GadgetMessageAcquisition>(mbq);
			GadgetContainerMessage< hoNDArray< std::complex<float> > >* daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq->cont());

			if (!acq || !daq) {
				GADGET_DEBUG1("Unable to interpret data on message Q\n");
				return GADGET_FAIL;
			}

			for (unsigned int c = 0; c < num_batches; c++) {
				float_complext* data_ptr = data_host->get_data_ptr();
				data_ptr += c*samples_per_interleave_*profiles_buffered+p*samples_per_interleave_;

				std::complex<float>* r_ptr = daq->getObjectPtr()->get_data_ptr();
				r_ptr += c*daq->getObjectPtr()->get_size(0);

				memcpy(data_ptr,r_ptr,samples_per_interleave_*sizeof(float_complext));
			}

			floatd2* traj_ptr = traj_host->get_data_ptr();
			traj_ptr += p*samples_per_interleave_;

			floatd2* t_ptr = host_traj_->get_data_ptr();
			t_ptr += acq->getObjectPtr()->idx.line*samples_per_interleave_;

			memcpy(traj_ptr,t_ptr,samples_per_interleave_*sizeof(floatd2));

			float* dcw_ptr = dcw_host->get_data_ptr();
			dcw_ptr += p*samples_per_interleave_;

			float* d_ptr = host_weights_->get_data_ptr();
			d_ptr += acq->getObjectPtr()->idx.line*samples_per_interleave_;

			memcpy(dcw_ptr,d_ptr,samples_per_interleave_*sizeof(float));

			mbq->release();
		}


		if (buffer_.message_count()) {
			GADGET_DEBUG1("Error occured, all messages should have been cleared off the buffer by now.\n");
			return GADGET_FAIL;
		}

		GadgetContainerMessage< CGSenseJob >* m4 =
				new GadgetContainerMessage< CGSenseJob >();

		m4->getObjectPtr()->dat_host_ = data_host;
		m4->getObjectPtr()->csm_host_ = csm_host;
		m4->getObjectPtr()->reg_host_ = reg_host;
		m4->getObjectPtr()->tra_host_ = traj_host;
		m4->getObjectPtr()->dcw_host_ = dcw_host;

		/*
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m4 =
				new GadgetContainerMessage< hoNDArray< std::complex<float> > >();


		if (!m4->getObjectPtr()->create(&image_dimensions_)) {
			GADGET_DEBUG1("Unable to allocate memory for combined image\n");
			m4->release();
			return GADGET_FAIL;
		}

		unsigned int npixels = image_dimensions_[0]*image_dimensions_[1];
		//std::complex<float>* recon_ptr    = reinterpret_cast< std::complex<float>* >(image_host->get_data_ptr());
		std::complex<float>* comb_ptr     = reinterpret_cast< std::complex<float>* >(m4->getObjectPtr()->get_data_ptr());
		memcpy(comb_ptr,reg_host.get()->get_data_ptr(),npixels*sizeof(float)*2);
		 */

		GadgetContainerMessage<GadgetMessageImage>* m3 =
				new GadgetContainerMessage<GadgetMessageImage>();

		m3->cont(m4);

		m3->getObjectPtr()->matrix_size[0] = image_dimensions_[0];
		m3->getObjectPtr()->matrix_size[1] = image_dimensions_[1];
		m3->getObjectPtr()->matrix_size[2] = 1;
		m3->getObjectPtr()->channels       = 1;
		m3->getObjectPtr()->data_idx_min       = m1->getObjectPtr()->min_idx;
		m3->getObjectPtr()->data_idx_max       = m1->getObjectPtr()->max_idx;
		m3->getObjectPtr()->data_idx_current   = m1->getObjectPtr()->idx;

		memcpy(m3->getObjectPtr()->position,m1->getObjectPtr()->position,
				sizeof(float)*3);

		memcpy(m3->getObjectPtr()->quarternion,m1->getObjectPtr()->quarternion,
				sizeof(float)*4);

		m3->getObjectPtr()->table_position = m1->getObjectPtr()->table_position;

		m3->getObjectPtr()->image_format = GADGET_IMAGE_COMPLEX_FLOAT; 
		m3->getObjectPtr()->image_index = ++image_counter_; 
		m3->getObjectPtr()->image_series_index = image_series_;

		if (this->next()->putq(m3) < 0) {
			m3->release();
			return GADGET_FAIL;
		}

	}
	//m1->release();
	return GADGET_OK;
}


GADGET_FACTORY_DECLARE(SpiralGadgetSW)
