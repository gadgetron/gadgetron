#include "NoiseAdjustGadget.h"
#include "Gadgetron.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "hoArmadillo.h"
#include "hoNDArray_elemwise.h"

namespace Gadgetron{

  NoiseAdjustGadget::NoiseAdjustGadget()
  : noise_decorrelation_calculated_(false)
  , number_of_noise_samples_(0)
  , noise_bw_scale_factor_(1.0f)
  , noise_dwell_time_us_(0.0f)
  , is_configured_(false)
  {
  }

  int NoiseAdjustGadget::process_config(ACE_Message_Block* mb)
  {
    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));
 
    receiver_noise_bandwidth_ = cfg->acquisitionSystemInformation().get().relativeReceiverNoiseBandwidth().present() ?
      cfg->acquisitionSystemInformation().get().relativeReceiverNoiseBandwidth().get() : 1.0;

    return GADGET_OK;
  }

  int NoiseAdjustGadget
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
	std::vector<unsigned long long> dims(2, channels);
			
	try{noise_covariance_matrix_.create(&dims);}
	catch (std::runtime_error& err)	{
	  GADGET_DEBUG_EXCEPTION(err, "Unable to allocate storage for noise covariance matrix\n" );
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
	  noise_bw_scale_factor_ = std::sqrt(2*acquisition_dwell_time_us_/noise_dwell_time_us_*receiver_noise_bandwidth_);
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
	  
	  std::vector<unsigned long long> dims(2, channels);
	  try{noise_covariance_matrixf_.create(&dims);}
	  catch (std::runtime_error& err){
	    GADGET_DEBUG_EXCEPTION(err,"Unable to allocate storage for noise covariance matrix (float)\n");
	    return GADGET_FAIL;
	  }
	  
	  // Armadillo can best do its template magic when we concatenate all the operations...
	  // 1. scale for number of samples
	  // 2. Cholesky decomposition
	  // 3. Invert lower triangular
	  // 4. Scale for noise BW

	  arma::cx_mat noise_cov = as_arma_matrix(&noise_covariance_matrix_);	  
	  arma::cx_fmat noise_covf = as_arma_matrix(&noise_covariance_matrixf_);

	  {	  
	    noise_covf = arma::conv_to<arma::cx_fmat>::from
	      (noise_bw_scale_factor_*arma::inv(arma::trimatu(arma::chol(noise_cov/number_of_noise_samples_))));
	  }
	  
	  noise_decorrelation_calculated_ = true;
	}
		
	if (noise_decorrelation_calculated_) {
	  arma::cx_fmat noise_covf = as_arma_matrix(&noise_covariance_matrixf_);
	  arma::cx_fmat am2 = as_arma_matrix(m2->getObjectPtr());	  
	  am2 = am2*arma::trimatu(noise_covf);
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
  
  GADGET_FACTORY_DECLARE(NoiseAdjustGadget)
  
} // namespace Gadgetron
