#include "MRINoiseAdjustGadget.h"
#include "Gadgetron.h"

#include "hoNDArray_fileio.h"
#include "matrix_vector_op.h"
#include "matrix_decomposition.h"
#include "GadgetronTimer.h"

#include "GadgetIsmrmrdReadWrite.h"

MRINoiseAdjustGadget::MRINoiseAdjustGadget()
: noise_decorrelation_calculated_(false)
, number_of_noise_samples_(0)
, noise_bw_scale_factor_(1.0f)
, is_configured_(false)
{
}


int MRINoiseAdjustGadget::process_config(ACE_Message_Block* mb)
{


	boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));


	receiver_noise_bandwidth_ = cfg->acquisitionSystemInformation().get().relativeReceiverNoiseBandwidth().present() ?
								cfg->acquisitionSystemInformation().get().relativeReceiverNoiseBandwidth().get() : 1.0;



	return GADGET_OK;
}

int MRINoiseAdjustGadget
::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

	bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
	unsigned int channels = m1->getObjectPtr()->active_channels;
	unsigned int samples = m1->getObjectPtr()->number_of_samples;

	if (is_noise) {
		noise_dwell_time_us_ = m1->getObjectPtr()->sample_time_us;

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


		/*
		for (unsigned int s = 0; s < samples; s++) {
			for (unsigned int i = 0; i < channels; i++) {
				for (unsigned int j = 0; j < channels; j++) {
					cc_ptr[j*channels + i] += (data_ptr[i * samples + s] * conj(data_ptr[j * samples + s]));
				}
			}
			number_of_noise_samples_++;
		}
		*/
		for (unsigned int s = 0; s < samples; s++) {
			for (unsigned int i = 0; i < channels; i++) {
				for (unsigned int j = 0; j < channels; j++) {
					cc_ptr[i*channels + j] += (data_ptr[i * samples + s] * conj(data_ptr[j * samples + s]));
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
				noise_bw_scale_factor_ = sqrt(2*acquisition_dwell_time_us_/noise_dwell_time_us_*receiver_noise_bandwidth_);
			}

			GADGET_DEBUG2("Noise dwell time: %f\n", noise_dwell_time_us_);
			GADGET_DEBUG2("Acquisition dwell time: %f\n", acquisition_dwell_time_us_);
			GADGET_DEBUG2("receiver_noise_bandwidth: %f\n", receiver_noise_bandwidth_);
			GADGET_DEBUG2("noise_bw_scale_factor: %f\n", noise_bw_scale_factor_);
			is_configured_ = true;
		}

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
				hoNDArray_choldc(&noise_covariance_matrix_);
				//choldc(cc_ptr, channels);

				//write_nd_array(&noise_covariance_matrix_, "CC_chol.cplx");

				//3. Invert lower triangular
				//inv_L(cc_ptr, channels);
				hoNDArray_inv_lower_triangular(&noise_covariance_matrix_);

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
				//static int data_written = 0;

				std::complex<float> alpha(1.0,0);
				if (hoNDArray_trmm(&noise_covariance_matrixf_, m2->getObjectPtr(), alpha) < 0) {
					GADGET_DEBUG1("Noise Decorrelation Failed\n");
					return GADGET_FAIL;
				}

				//Noise decorrelate
				/*
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
				/*
				/*
				if (!data_written) {
					write_nd_array(&noise_covariance_matrixf_, "noise_decorr_matrix.cplx");
					write_nd_array(m2->getObjectPtr(), "data_nodcx.cplx");
					write_nd_array(m3->getObjectPtr(), "data_dcx.cplx");
					data_written++;
				}
				*/

				//m1->cont(m3);
				//m2->release();

				/*
				if (!noise_decorrelation(m2->getObjectPtr()->get_data_ptr(), samples, channels, noise_covariance_matrix_.get_data_ptr())) {
					GADGET_DEBUG1("Noise Decorrelation Failed\n");
					return GADGET_FAIL;
				}
				if (!data_written) {
					write_nd_array(&noise_covariance_matrixf_, "noise_decorr_matrix.cplx");
					write_nd_array(m2->getObjectPtr(), "data_dcx.cplx");
					data_written++;
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
