#include "MRINoiseAdjustGadget.h"
#include "Gadgetron.h"
#include "GadgetXml.h"

#include "hoNDArray_fileio.h"
#include "matrix_vector_op.h"

void choldc(std::complex<double> *a, int n)
{
	int i,j,k;

	for (k= 0; k < n; k++)
	{
		a[k*n+k] = std::complex<double>(sqrt(real(a[k*n+k])),0.0);

		for (i = k+1; i < n; i++)
		{
			a[k*n+i] = a[k*n+i]/a[k*n+k];
		}

		for (j = k + 1; j < n; j++)
		{
			for (i = j; i < n; i++)
			{
				a[j*n+i] -= conj(a[k*n+j])*a[k*n+i];
			}
		}
	}
}

void inv_L(std::complex<double> *a, int n)
{
	int i,j,k;

	std::complex<double> sum;

	for (i = 0; i < n; i++)
	{

		a[i*n+i] = std::complex<double>(1.0/real(a[i*n+i]),0.0);
		for (j = i+1; j < n; j++)
		{
			sum = std::complex<double>(0.0,0.0);
			for (k = i; k < j; k++)
			{
				sum -= a[k*n+j]*a[i*n+k];
			}
			a[i*n+j] = sum/a[j*n+j];
		}
	}
}

bool noise_decorrelation(std::complex<float>* data, int elements, int coils, std::complex<double>* inv_L_psi)
{
	int i,j,k;

	/* We need some temporary storrage to store the data for one element before overwriting the original data */
	std::complex<double>* tmp_data = new std::complex<double>[coils];

	if (tmp_data == 0)
	{
		return false;
	}

	for (i = 0; i < elements; i++)
	{
		for (j = 0; j < coils; j++)
		{
			tmp_data[j] = std::complex<double>(0.0,0.0);
		}

		for (j = 0; j < coils; j++)
		{
			for (k = 0; k <= j; k++)
			{
				tmp_data[j] += inv_L_psi[k*coils+j] * static_cast< std::complex<double> >(data[k*elements+i]);
			}
		}

		for (j = 0; j < coils; j++)
		{
			data[j*elements+i] = tmp_data[j];
		}
	}

	/* Clean up */
	delete [] tmp_data;

	return true;
}



MRINoiseAdjustGadget::MRINoiseAdjustGadget()
: noise_decorrelation_calculated_(false)
, number_of_noise_samples_(0)
, noise_bw_scale_factor_(1.0f)
{

}


int MRINoiseAdjustGadget::process_config(ACE_Message_Block* mb)
{

	TiXmlDocument doc;
	doc.Parse(mb->rd_ptr());

	GadgetXMLNode n = GadgetXMLNode(&doc);

	noise_dwell_time_us_ = 7680.0f/(n.get<double>(std::string("gadgetron.encoding.kspace.readout_length.value"))[0]);//n.get<double>(std::string("noise_dwell_time_us.value"))[0]; //GetDoubleParameterValueFromXML(&doc, "encoding", "noise_dwell_time_us");
	acquisition_dwell_time_us_ = (n.get<double>(std::string("gadgetron.encoding.acquisition_dwell_time_ns.value"))[0])/1000.0f; //GetDoubleParameterValueFromXML(&doc, "encoding", "acquisition_dwell_time_us");
	receiver_noise_bandwidth_ = 0.79;//n.get<double>(std::string("receiver_noise_bandwidth.value"))[0];//GetDoubleParameterValueFromXML(&doc, "encoding", "receiver_noise_bandwidth");

	if ((noise_dwell_time_us_ == 0.0f) || (acquisition_dwell_time_us_ == 0.0f)) {
		noise_bw_scale_factor_ = 1.0f;
	} else {
		noise_bw_scale_factor_ = sqrt(2*acquisition_dwell_time_us_/noise_dwell_time_us_*receiver_noise_bandwidth_);
	}

	GADGET_DEBUG2("Noise dwell time: %f\n", noise_dwell_time_us_);
	GADGET_DEBUG2("Acquisition dwell time: %f\n", acquisition_dwell_time_us_);
	GADGET_DEBUG2("receiver_noise_bandwidth: %f\n", receiver_noise_bandwidth_);
	GADGET_DEBUG2("noise_bw_scale_factor: %f\n", noise_bw_scale_factor_);


	return GADGET_OK;
}

int MRINoiseAdjustGadget
::process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

	bool is_noise = m1->getObjectPtr()->flags & GADGET_FLAG_IS_NOISE_SCAN;
	unsigned int channels = m1->getObjectPtr()->channels;
	unsigned int samples = m1->getObjectPtr()->samples;

	if (is_noise) {
		//If noise covariance matrix is not allocated
		if (noise_covariance_matrix_.get_number_of_elements() != channels*channels) {
			std::vector<unsigned int> dims(2, channels);
			if (!noise_covariance_matrix_.create(&dims)) {
				GADGET_DEBUG1("Unable to allocate storage for noise covariance matrix\n");
				return GADGET_FAIL;
			} else {
				noise_covariance_matrix_.clear(std::complex<double>(0.0,0.0));
			}
			number_of_noise_samples_ = 0;
		}

		std::complex<double>* cc_ptr = noise_covariance_matrix_.get_data_ptr();
		std::complex<float>* data_ptr = m2->getObjectPtr()->get_data_ptr();


		for (unsigned int s = 0; s < samples; s++) {
			for (unsigned int i = 0; i < channels; i++) {
				for (unsigned int j = 0; j < channels; j++) {
					cc_ptr[j*channels + i] += (data_ptr[i * samples + s] * conj(data_ptr[j * samples + s]));
				}
			}
			number_of_noise_samples_++;
		}
	} else {
		if (number_of_noise_samples_ > 0) {
			if (!noise_decorrelation_calculated_) {
				GADGET_DEBUG1("Calculating noise decorrelation\n");
				//1. scale for number of samples
				std::complex<double>* cc_ptr = noise_covariance_matrix_.get_data_ptr();
				for (unsigned int i = 0; i < channels*channels; i++) {
					cc_ptr[i] /= number_of_noise_samples_;
				}

				//write_nd_array(&noise_covariance_matrix_, "CC.cplx");

				//2. Cholesky decomposition
				choldc(cc_ptr, channels);

				//write_nd_array(&noise_covariance_matrix_, "CC_chol.cplx");

				//3. Invert lower triangular
				inv_L(cc_ptr, channels);

				//write_nd_array(&noise_covariance_matrix_, "CC_chol_inv_L.cplx");

				//4. Scale for noise BW
				for (unsigned int i = 0; i < channels*channels; i++) {
					cc_ptr[i] *= noise_bw_scale_factor_;
				}

				/* Copy to float precision */
				std::vector<unsigned int> dims(2, channels);
				if (!noise_covariance_matrixf_.create(&dims)) {
					GADGET_DEBUG1("Unable to allocate storage for noise covariance matrix (float)\n");
					return GADGET_FAIL;
				} else {
					noise_covariance_matrixf_.clear(std::complex<float>(0.0,0.0));
				}

				std::complex<float>* ccf_ptr = noise_covariance_matrixf_.get_data_ptr();
				for (unsigned int i = 0; i < channels*channels; i++) {
					ccf_ptr[i] = cc_ptr[i];
				}

				noise_decorrelation_calculated_ = true;
			}

			if (noise_decorrelation_calculated_) {
				//Noise decorrelate
				GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 =
						new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

				if (!m3->getObjectPtr()->create(m2->getObjectPtr()->get_dimensions().get())) {
					GADGET_DEBUG1("Unable to allocate storage for decorrelated data\n");
					return GADGET_FAIL;
				}
				std::complex<float> alpha(1.0,0);
				std::complex<float> beta(0.0,0);
				if (hoNDArray_gemm(&noise_covariance_matrixf_, m2->getObjectPtr(), alpha, m3->getObjectPtr(), beta) < 0) {
					GADGET_DEBUG1("Noise Decorrelation Failed\n");
					return GADGET_FAIL;
				}
				m1->cont(m3);
				m2->release();
				/*
				if (!noise_decorrelation(m2->getObjectPtr()->get_data_ptr(), samples, channels, noise_covariance_matrix_.get_data_ptr())) {
					GADGET_DEBUG1("Noise Decorrelation Failed\n");
					return GADGET_FAIL;
				}
				*/
			}
		}
		//It is enough to put the first one, since they are linked
		if (this->next()->putq(m1) == -1) {
			ACE_ERROR_RETURN( (LM_ERROR,
					ACE_TEXT("%p\n"),
					ACE_TEXT("NoiseAdjustGadget::process, passing data on to next gadget")),
					-1);
		}

	}

	return GADGET_OK;
}


GADGET_FACTORY_DECLARE(MRINoiseAdjustGadget)
