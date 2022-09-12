#include "NoiseAdjustGadget_unoptimized.h"
#include "hoNDArray_fileio.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{

void choldc(std::complex<double> *a, int n)
{
	int i,j,k;

	for (k= 0; k < n; k++)
	{
		a[k*n+k] = std::complex<double>(std::sqrt(real(a[k*n+k])),0.0);

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



NoiseAdjustGadget_unoptimized::NoiseAdjustGadget_unoptimized()
: noise_decorrelation_calculated_(false)
, number_of_noise_samples_(0)
, noise_bw_scale_factor_(1.0f)
, is_configured_(false)
{

}


int NoiseAdjustGadget_unoptimized::process_config(ACE_Message_Block* mb)
{
  ISMRMRD::IsmrmrdHeader h;
  ISMRMRD::deserialize(mb->rd_ptr(),h);
  
  if ( h.acquisitionSystemInformation ) {
    receiver_noise_bandwidth_ = (float)(h.acquisitionSystemInformation->relativeReceiverNoiseBandwidth ?
					*h.acquisitionSystemInformation->relativeReceiverNoiseBandwidth : 1.0f);
    
    GDEBUG_STREAM("receiver_noise_bandwidth_ is " << receiver_noise_bandwidth_);
  }
  
  return GADGET_OK;
}

int NoiseAdjustGadget_unoptimized
::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

	bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
	unsigned int channels = m1->getObjectPtr()->active_channels;
	unsigned int samples = m1->getObjectPtr()->number_of_samples;

	if (is_noise) {
		noise_dwell_time_us_ = m1->getObjectPtr()->sample_time_us;
		//If noise covariance matrix is not allocated
		if (noise_covariance_matrix_.get_number_of_elements() != channels*channels) {
			std::vector<size_t> dims(2, channels);
			try{ noise_covariance_matrix_.create(dims);}
			catch (std::runtime_error &err){
				GEXCEPTION(err,"Unable to allocate storage for noise covariance matrix\n");
				return GADGET_FAIL;
			}
			noise_covariance_matrix_.fill(std::complex<double>(0.0,0.0));

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
		acquisition_dwell_time_us_ = m1->getObjectPtr()->sample_time_us;
		if (!is_configured_) {
			if ((noise_dwell_time_us_ == 0.0f) || (acquisition_dwell_time_us_ == 0.0f)) {
				noise_bw_scale_factor_ = 1.0f;
			} else {
				noise_bw_scale_factor_ = std::sqrt(2*acquisition_dwell_time_us_/noise_dwell_time_us_*receiver_noise_bandwidth_);
			}

			GDEBUG("Noise dwell time: %f\n", noise_dwell_time_us_);
			GDEBUG("Acquisition dwell time: %f\n", acquisition_dwell_time_us_);
			GDEBUG("receiver_noise_bandwidth: %f\n", receiver_noise_bandwidth_);
			GDEBUG("noise_bw_scale_factor: %f\n", noise_bw_scale_factor_);
			is_configured_ = true;
		}
		if (number_of_noise_samples_ > 0) {
			if (!noise_decorrelation_calculated_) {
				GDEBUG("Calculating noise decorrelation\n");
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

				noise_decorrelation_calculated_ = true;
			}

			if (noise_decorrelation_calculated_) {
				//Noise decorrelate
				if (!noise_decorrelation(m2->getObjectPtr()->get_data_ptr(), samples, channels, noise_covariance_matrix_.get_data_ptr())) {
					GDEBUG("Noise Decorrelation Failed\n");
					return GADGET_FAIL;
				}
			}
		}
		//It is enough to put the first one, since they are linked
		if (this->next()->putq(m1) == -1) {
		  GERROR("NoiseAdjustGadget_unoptimized::process, passing data on to next gadget");
		  return -1;
		}

	}

	return GADGET_OK;
}


GADGET_FACTORY_DECLARE(NoiseAdjustGadget_unoptimized)
}
