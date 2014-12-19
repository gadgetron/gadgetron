#include "NoiseAdjustGadget.h"
#include "Gadgetron.h"
#include "hoArmadillo.h"
#include "hoNDArray_elemwise.h"
#include "GadgetronCommon.h"
#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_elemwise.h"

#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

#ifndef _WIN32
#include <sys/types.h>
#include <sys/stat.h>
#endif // _WIN32

namespace Gadgetron{

  NoiseAdjustGadget::NoiseAdjustGadget()
    : noise_decorrelation_calculated_(false)
    , number_of_noise_samples_(0)
    , number_of_noise_samples_per_acquisition_(0)
    , noise_bw_scale_factor_(1.0f)
    , noise_dwell_time_us_(-1.0f)
    , is_configured_(false)
    , computed_in_close_(false)
    , noiseCovarianceLoaded_(false)
  {
    noise_dependency_prefix_ = "GadgetronNoiseCovarianceMatrix";

    measurement_id_.clear();
    measurement_id_of_noise_dependency_.clear();

    noise_dwell_time_us_preset_ = 5;

    perform_noise_adjust_ = true;

    gt_timer_.set_timing_in_destruction(false);
    performTiming_ = false;
  }

  NoiseAdjustGadget::~NoiseAdjustGadget()
  {

  }

  int NoiseAdjustGadget::process_config(ACE_Message_Block* mb)
  {

    boost::shared_ptr<std::string> str = this->get_string_value("workingDirectory");

    if ( !str->empty() ) {
      noise_dependency_folder_ = *str;
    }
    else {
#ifdef _WIN32
      noise_dependency_folder_ = std::string("c:\\temp\\gadgetron\\");
#else
      noise_dependency_folder_ =  std::string("/tmp/gadgetron/");
#endif // _WIN32
    }

    GADGET_DEBUG2("Folder to store noise dependencies is %s\n", noise_dependency_folder_.c_str());

    str = this->get_string_value("noise_dependency_prefix");
    if ( !str->empty() ) noise_dependency_prefix_ = *str;

    performTiming_ = this->get_bool_value("performTiming");

    perform_noise_adjust_ = this->get_string_value("perform_noise_adjust")->size() ? this->get_bool_value("perform_noise_adjust") : true;

    GADGET_DEBUG2("NoiseAdjustGadget::perform_noise_adjust_ is %d\n", perform_noise_adjust_);

    noise_dwell_time_us_preset_ = (float)this->get_double_value("noise_dwell_time_us_preset");
    if ( noise_dwell_time_us_preset_ == 0 ) noise_dwell_time_us_preset_ = 5;

    ISMRMRD::deserialize(mb->rd_ptr(),current_ismrmrd_header_);
    
    if ( current_ismrmrd_header_.acquisitionSystemInformation )
      {
	receiver_noise_bandwidth_ = (float)(current_ismrmrd_header_.acquisitionSystemInformation->relativeReceiverNoiseBandwidth ?
					    *current_ismrmrd_header_.acquisitionSystemInformation->relativeReceiverNoiseBandwidth : 0.793f);

	GADGET_MSG("receiver_noise_bandwidth_ is " << receiver_noise_bandwidth_);
      }

    // find the measurementID of this scan
    if ( current_ismrmrd_header_.measurementInformation )
      {
	if ( current_ismrmrd_header_.measurementInformation->measurementID )
	  {
	    measurement_id_ = *current_ismrmrd_header_.measurementInformation->measurementID;
	    GADGET_MSG("Measurement ID is " << measurement_id_);
	  }

	// find the noise depencies if any
	if ( current_ismrmrd_header_.measurementInformation->measurementDependency.size() > 0 )
	  {
	    measurement_id_of_noise_dependency_.clear();

	    std::vector<ISMRMRD::MeasurementDependency>::const_iterator iter = current_ismrmrd_header_.measurementInformation->measurementDependency.begin();
	    for ( ; iter!= current_ismrmrd_header_.measurementInformation->measurementDependency.end(); iter++ )
	      {
		std::string dependencyType = iter->dependencyType;
		std::string dependencyID = iter->measurementID;

		GADGET_MSG("Found dependency measurement : " << dependencyType << " with ID " << dependencyID);
            
		if ( dependencyType=="Noise" || dependencyType=="noise" ) {
		  measurement_id_of_noise_dependency_ = dependencyID;
		}
	      }
        
	    if ( !measurement_id_of_noise_dependency_.empty() ) {
	      GADGET_MSG("Measurement ID of noise dependency is " << measurement_id_of_noise_dependency_);
		  
	      full_name_stored_noise_dependency_ = this->generateNoiseDependencyFilename(generateMeasurementIdOfNoiseDependency(measurement_id_of_noise_dependency_));
	      GADGET_MSG("Stored noise dependency is " << full_name_stored_noise_dependency_);
		  
	      // try to load the precomputed noise prewhitener
	      if ( !this->loadNoiseCovariance() )
		{
		  GADGET_MSG("Stored noise dependency is NOT found : " << full_name_stored_noise_dependency_);
		  noiseCovarianceLoaded_ = false;
		  noise_dwell_time_us_ = -1;
		  noise_covariance_matrixf_.clear();
		}
	      else
		{
		  GADGET_MSG("Stored noise dependency is found : " << full_name_stored_noise_dependency_);
		  GADGET_MSG("Stored noise dwell time in us is " << noise_dwell_time_us_);
		  GADGET_MSG("Stored noise channel number is " << noise_covariance_matrixf_.get_size(0));
		  noiseCovarianceLoaded_ = true;
		  number_of_noise_samples_ = 1; //When we load the matrix, it is already scaled.
		}
	    }
	  }
      }

    // limit the number of threads used to be 1
#ifdef USE_OMP
    omp_set_num_threads(1);
    GADGET_MSG("NoiseAdjustGadget:omp_set_num_threads(1) ... ");
#endif // USE_OMP

    return GADGET_OK;
  }

  std::string NoiseAdjustGadget::generateMeasurementIdOfNoiseDependency(const std::string& noise_id)
  {
    // find the scan prefix
    std::string measurementStr = measurement_id_;
    size_t ind  = measurement_id_.find_last_of ("_");
    if ( ind != std::string::npos ) {
      measurementStr = measurement_id_.substr(0, ind);
      measurementStr.append("_");
      measurementStr.append(noise_id);
    }
   
    return measurementStr;
  }

  std::string NoiseAdjustGadget::generateNoiseDependencyFilename(const std::string& measurement_id)
  {
    std::string full_name_stored_noise_dependency;

    full_name_stored_noise_dependency = noise_dependency_folder_;
    full_name_stored_noise_dependency.append("/");
    full_name_stored_noise_dependency.append(noise_dependency_prefix_);
    full_name_stored_noise_dependency.append("_");
    full_name_stored_noise_dependency.append(measurement_id);

    return full_name_stored_noise_dependency;
  }

  bool NoiseAdjustGadget::loadNoiseCovariance()
  {
    std::ifstream infile;
    infile.open (full_name_stored_noise_dependency_.c_str(), std::ios::in|std::ios::binary);

    if (infile.good() )
      {
	//Read the XML header of the noise scan
	uint32_t xml_length;
	infile.read( reinterpret_cast<char*>(&xml_length), 4);
	std::string xml_str(xml_length,'\0');
	infile.read(const_cast<char*>(xml_str.c_str()), xml_length);
	ISMRMRD::deserialize(xml_str.c_str(), noise_ismrmrd_header_);
	
	infile.read( reinterpret_cast<char*>(&noise_dwell_time_us_), sizeof(float));

	size_t len;
	infile.read( reinterpret_cast<char*>(&len), sizeof(size_t));

	char* buf = new char[len];
	if ( buf == NULL ) return false;

	infile.read(buf, len);

	if ( !noise_covariance_matrixf_.deserialize(buf, len) )
	  {
	    delete [] buf;
	    return false;
	  }

	delete [] buf;
	infile.close();
      }
    else
      {
	GADGET_ERROR_MSG("Noise prewhitener file is not good for reading");
	return false;
      }

    return true;
  }

  bool NoiseAdjustGadget::saveNoiseCovariance()
  {
    char* buf = NULL;
    size_t len(0);

    //Scale the covariance matrix before saving
    hoNDArray< std::complex<float> > covf(noise_covariance_matrixf_);

    if (number_of_noise_samples_ > 1) {
      covf *= std::complex<float>(1.0/(float)(number_of_noise_samples_-1),0.0);
    }

    if ( !covf.serialize(buf, len) )
      {
	GADGET_ERROR_MSG("Noise prewhitener serialization failed ...");
	return false;
      }

    std::stringstream xml_ss;
    ISMRMRD::serialize(current_ismrmrd_header_, xml_ss);
    std::string xml_str = xml_ss.str();
    uint32_t xml_length = static_cast<uint32_t>(xml_str.size());

    std::ofstream outfile;
    std::string filename  = this->generateNoiseDependencyFilename(measurement_id_);
    outfile.open (filename.c_str(), std::ios::out|std::ios::binary);

    if (outfile.good())
      {
	GADGET_MSG("write out the noise dependency file : " << filename);
	outfile.write( reinterpret_cast<char*>(&xml_length), 4);
	outfile.write( xml_str.c_str(), xml_length );
	outfile.write( reinterpret_cast<char*>(&noise_dwell_time_us_), sizeof(float));
	outfile.write( reinterpret_cast<char*>(&len), sizeof(size_t));
	outfile.write(buf, len);
	outfile.close();

	// set the permission for the noise file to be rewritable
#ifndef _WIN32
	int res = chmod(filename.c_str(), S_IRUSR|S_IWUSR|S_IXUSR|S_IRGRP|S_IWGRP|S_IXGRP|S_IROTH|S_IWOTH|S_IXOTH);
	if ( res != 0 )
	  {
	    GADGET_ERROR_MSG("Changing noise prewhitener file permission failed ...");
	  }
#endif // _WIN32
      }
    else
      {
	delete [] buf;
	GADGET_ERROR_MSG("Noise prewhitener file is not good for writing");
	return false;
      }

    delete [] buf;
    return true;
  }

  void NoiseAdjustGadget::computeNoisePrewhitener()
  {
    GADGET_START_TIMING_CONDITION(gt_timer_, "compute noise prewhitener ... ", performTiming_);

    if ( noise_dwell_time_us_ > 0 )  {
      GADGET_MSG("Noise dwell time: " << noise_dwell_time_us_);
      GADGET_MSG("receiver_noise_bandwidth: " << receiver_noise_bandwidth_);
      
      if (!noise_decorrelation_calculated_) {
	GADGET_MSG("Calculating noise decorrelation");
	
	// Armadillo can best do its template magic when we concatenate all the operations...
	// 1. scale for number of samples
	// 2. Cholesky decomposition
	// 3. Invert lower triangular
	
	noise_prewhitener_matrixf_ = noise_covariance_matrixf_;
	arma::cx_fmat noise_covf = as_arma_matrix(&noise_prewhitener_matrixf_);
	
	noise_covf = arma::inv(arma::trimatu(arma::chol(noise_covf)));
	noise_decorrelation_calculated_ = true;
      }
    }
      
    GADGET_STOP_TIMING_CONDITION(gt_timer_, performTiming_);
  }

  int NoiseAdjustGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
    unsigned int channels = m1->getObjectPtr()->active_channels;
    unsigned int samples = m1->getObjectPtr()->number_of_samples;

    //TODO: Remove this
    if ( measurement_id_.empty() ) {
      unsigned int muid = m1->getObjectPtr()->measurement_uid;
      std::ostringstream ostr;
      ostr << muid;
      measurement_id_ = ostr.str();
    }

    if ( is_noise ) {
      // this noise can be from a noise scan or it can be from the built-in noise
      if ( number_of_noise_samples_per_acquisition_ == 0 ) {
	number_of_noise_samples_per_acquisition_ = samples;
      }

      //TODO: This should probably just be:
      //if (noise_dwell_time_us_ < 0) noise_dwell_time_us_ = m1->getObjectPtr()->sample_time_us;
      if ( noise_dwell_time_us_ < 0 ) {
	//TODO: Investigate this surface coil correction business. 
	//if ( !is_scc_correction && number_of_noise_samples_per_acquisition_>0 ) { IS THIS NECESSARY??
	if ( number_of_noise_samples_per_acquisition_>0 ) {
	  noise_dwell_time_us_ = m1->getObjectPtr()->sample_time_us;
	} else {
	  noise_dwell_time_us_ = noise_dwell_time_us_preset_;
	}
      }

      //If noise covariance matrix is not allocated
      if (noise_covariance_matrixf_.get_number_of_elements() != channels*channels) {
	std::vector<size_t> dims(2, channels);
	try {
	  noise_covariance_matrixf_.create(&dims);
	  noise_covariance_matrixf_once_.create(&dims);
	} catch (std::runtime_error& err) {
	  GADGET_DEBUG_EXCEPTION(err, "Unable to allocate storage for noise covariance matrix\n" );
	  return GADGET_FAIL;
	}

	Gadgetron::clear(noise_covariance_matrixf_);
	Gadgetron::clear(noise_covariance_matrixf_once_);
	number_of_noise_samples_ = 0;
      }

      std::complex<float>* cc_ptr = noise_covariance_matrixf_.get_data_ptr();
      std::complex<float>* data_ptr = m2->getObjectPtr()->get_data_ptr();
      
      readout_ = *m2->getObjectPtr();
      gemm(noise_covariance_matrixf_once_, readout_, true, *m2->getObjectPtr(), false);
      Gadgetron::add(noise_covariance_matrixf_once_, noise_covariance_matrixf_, noise_covariance_matrixf_);
      
      number_of_noise_samples_ += samples;
      m1->release();
      return GADGET_OK;
    }


    //We should only reach this code if this data is not noise.
    if ( perform_noise_adjust_ ) {
      //Calculate the prewhitener if it has not been done
      if (!noise_decorrelation_calculated_) {
	if (number_of_noise_samples_ > 1) {
	  //Scale
	  noise_covariance_matrixf_ *= std::complex<float>(1.0/(float)(number_of_noise_samples_-1));
	  number_of_noise_samples_ = 1; //Scaling has been done
	}
	computeNoisePrewhitener();
	acquisition_dwell_time_us_ = m1->getObjectPtr()->sample_time_us;
	if ((noise_dwell_time_us_ == 0.0f) || (acquisition_dwell_time_us_ == 0.0f)) {
	  noise_bw_scale_factor_ = 1.0f;
	}
	else {
	  noise_bw_scale_factor_ = (float)std::sqrt(2.0*acquisition_dwell_time_us_/noise_dwell_time_us_*receiver_noise_bandwidth_);
	}

	noise_prewhitener_matrixf_ *= std::complex<float>(noise_bw_scale_factor_,0.0);

	GADGET_MSG("Noise dwell time: " << noise_dwell_time_us_);
	GADGET_MSG("Acquisition dwell time:" << acquisition_dwell_time_us_);
	GADGET_MSG("receiver_noise_bandwidth: " << receiver_noise_bandwidth_);
	GADGET_MSG("noise_bw_scale_factor: " << noise_bw_scale_factor_);
      }

      
      //Apply prewhitener
      hoNDArray<std::complex<float> > tmp(*m2->getObjectPtr());
      gemm(*m2->getObjectPtr(), tmp, noise_prewhitener_matrixf_);
      
    }

    if (this->next()->putq(m1) == -1) {
      GADGET_DEBUG1("Error passing on data to next gadget\n");
      return GADGET_FAIL;
    }
    
    return GADGET_OK;

  }
      /*
      {
	if ( !is_noise )
	  {
	  }
	else
	  {
	    m1->release();
	  }

	return GADGET_OK;
      }



    if ( use_stored_noise_prewhitener_ )
      {
	if ( !is_noise )
	  {
	    acquisition_dwell_time_us_ = m1->getObjectPtr()->sample_time_us;
	    if (!is_configured_)
	      {
		if ((noise_dwell_time_us_ == 0.0f) || (acquisition_dwell_time_us_ == 0.0f))
		  {
		    noise_bw_scale_factor_ = 1.0f;
		  }
		else
		  {
		    noise_bw_scale_factor_ = (float)std::sqrt(2.0*acquisition_dwell_time_us_/noise_dwell_time_us_*receiver_noise_bandwidth_);
		  }

		GADGET_MSG("Noise dwell time: " << noise_dwell_time_us_);
		GADGET_MSG("Acquisition dwell time:" << acquisition_dwell_time_us_);
		GADGET_MSG("receiver_noise_bandwidth: " << receiver_noise_bandwidth_);
		GADGET_MSG("noise_bw_scale_factor: " << noise_bw_scale_factor_);
		is_configured_ = true;
	      }

	    if ( !noise_decorrelation_calculated_ )
	      {
		this->computeNoisePrewhitener();
		arma::cx_fmat noise_covf = as_arma_matrix(&noise_prewhitener_matrixf_);
		noise_covf *= noise_bw_scale_factor_;
		noise_decorrelation_calculated_ = true;
	      }

	    if (noise_decorrelation_calculated_)
	      {
		// GADGET_START_TIMING_CONDITION(gt_timer_, "apply noise prewhitener ... ", performTiming_);

		if ( data_prewhitened_.get_size(0)!=m2->getObjectPtr()->get_size(0) 
		     || data_prewhitened_.get_size(1)!=m2->getObjectPtr()->get_size(1) )
		  {
		    data_prewhitened_.create(m2->getObjectPtr()->get_dimensions());
		  }

		memcpy(data_prewhitened_.begin(), m2->getObjectPtr()->begin(), m2->getObjectPtr()->get_number_of_bytes());

		gemm(*m2->getObjectPtr(), data_prewhitened_, noise_covariance_matrixf_);

		// GADGET_STOP_TIMING_CONDITION(gt_timer_, performTiming_);
	      }

	    if (this->next()->putq(m1) == -1)
	      {
		ACE_ERROR_RETURN( (LM_ERROR,
				   ACE_TEXT("%p\n"),
				   ACE_TEXT("NoiseAdjustGadget::process, passing data on to next gadget")),
				  -1);
	      }
	  }
	else
	  {
	    m1->release();
	  }
      }
    else
      {
	if ( is_noise )
	  {
	    // this noise can be from a noise scan or it can be from the built-in noise
	    if ( number_of_noise_samples_per_acquisition_ == 0 )
	      {
		number_of_noise_samples_per_acquisition_ = samples;
	      }

	    if ( noise_dwell_time_us_ < 0 )
	      {
		if ( !is_scc_correction && number_of_noise_samples_per_acquisition_>0 )
		  {
		    noise_dwell_time_us_ = m1->getObjectPtr()->sample_time_us;
		  }
		else
		  {
		    noise_dwell_time_us_ = noise_dwell_time_us_preset_;
		  }
	      }

	    //If noise covariance matrix is not allocated
	    if (noise_covariance_matrixf_.get_number_of_elements() != channels*channels)
	      {
		std::vector<size_t> dims(2, channels);

		try
		  {
		    noise_covariance_matrixf_.create(&dims);
		    noise_covariance_matrixf_once_.create(&dims);
		  }
		catch (std::runtime_error& err)
		  {
		    GADGET_DEBUG_EXCEPTION(err, "Unable to allocate storage for noise covariance matrix\n" );
		    return GADGET_FAIL;
		  }
		Gadgetron::clear(noise_covariance_matrixf_);
		Gadgetron::clear(noise_covariance_matrixf_once_);
		number_of_noise_samples_ = 0;
	      }

	    std::complex<float>* cc_ptr = noise_covariance_matrixf_.get_data_ptr();
	    std::complex<float>* data_ptr = m2->getObjectPtr()->get_data_ptr();

	    readout_ = *m2->getObjectPtr();
	    gemm(noise_covariance_matrixf_once_, readout_, true, *m2->getObjectPtr(), false);
	    Gadgetron::add(noise_covariance_matrixf_once_, noise_covariance_matrixf_, noise_covariance_matrixf_);

	    number_of_noise_samples_ += samples;
	    m1->release();
	  }
	else
	  {
	    if ( noise_dwell_time_us_ > 0 )
	      {
		acquisition_dwell_time_us_ = m1->getObjectPtr()->sample_time_us;
		if (!is_configured_)
		  {
		    if ((noise_dwell_time_us_ == 0.0f) || (acquisition_dwell_time_us_ == 0.0f))
		      {
			noise_bw_scale_factor_ = 1.0f;
		      }
		    else
		      {
			noise_bw_scale_factor_ = (float)std::sqrt(2.0*acquisition_dwell_time_us_/noise_dwell_time_us_*receiver_noise_bandwidth_);
		      }

		    GADGET_MSG("Noise dwell time: " << noise_dwell_time_us_);
		    GADGET_MSG("Acquisition dwell time:" << acquisition_dwell_time_us_);
		    GADGET_MSG("receiver_noise_bandwidth: " << receiver_noise_bandwidth_);
		    GADGET_MSG("noise_bw_scale_factor: " << noise_bw_scale_factor_);
		    is_configured_ = true;
		  }

		if (number_of_noise_samples_ > 0)
		  {
		    if (!noise_decorrelation_calculated_)
		      {
			if ( is_scc_correction )
			  {
			    this->computeNoisePrewhitener();
			  }
			else
			  {
			    this->computeNoisePrewhitener();
			  }

			// apply the scaling
			arma::cx_fmat noise_covf = as_arma_matrix(&noise_covariance_matrixf_);
			noise_covf *= noise_bw_scale_factor_;
		      }
		    else
		      {
			if ( m2->getObjectPtr()->get_size(1) == noise_covariance_matrixf_.get_size(0) )
			  {
			    if ( data_prewhitened_.get_size(0)!=m2->getObjectPtr()->get_size(0) 
				 || data_prewhitened_.get_size(1)!=m2->getObjectPtr()->get_size(1) )
			      {
				data_prewhitened_.create(m2->getObjectPtr()->get_dimensions());
			      }

			    memcpy(data_prewhitened_.begin(), m2->getObjectPtr()->begin(), m2->getObjectPtr()->get_number_of_bytes());

			    gemm(*m2->getObjectPtr(), data_prewhitened_, noise_covariance_matrixf_);
			  }
		      }
		  }
	      }

	    //It is enough to put the first one, since they are linked
	    if (this->next()->putq(m1) == -1)
	      {
		ACE_ERROR_RETURN( (LM_ERROR,
				   ACE_TEXT("%p\n"),
				   ACE_TEXT("NoiseAdjustGadget::process, passing data on to next gadget")),
				  -1);
	      }
	  }
      }

    // GADGET_STOP_TIMING_CONDITION(gt_timer_, performTiming_);

    return GADGET_OK;
    }	  
      */

  int NoiseAdjustGadget::close(unsigned long flags)
  {
    if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

    static bool saved = false;

    if ( !noiseCovarianceLoaded_  && !saved ){
      if ( noise_dwell_time_us_ < 0 ) noise_dwell_time_us_ = noise_dwell_time_us_preset_; // this scan is a noise measurement
      saveNoiseCovariance();
      saved = true;
    }  

    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(NoiseAdjustGadget)

} // namespace Gadgetron
