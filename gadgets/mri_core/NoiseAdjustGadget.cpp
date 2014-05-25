#include "NoiseAdjustGadget.h"
#include "Gadgetron.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "hoArmadillo.h"
#include "hoNDArray_elemwise.h"
#include "GadgetronCommon.h"
#include "hoMatrix.h"
#include "hoMatrix_util.h"
#include "hoNDArray_math_util.h"

namespace Gadgetron{

    NoiseAdjustGadget::NoiseAdjustGadget()
        : noise_decorrelation_calculated_(false)
        , number_of_noise_samples_(0)
        , number_of_noise_samples_per_acquisition_(0)
        , noise_bw_scale_factor_(1.0f)
        , noise_dwell_time_us_(-1.0f)
        , is_configured_(false)
        , computed_in_close_(false)
        , use_stored_noise_prewhitener_(false)
    {
        noise_dependency_prefix_ = "GadgetronNoisePreWhitener";

        patient_id_.clear();
        study_id_.clear();
        measurement_id_.clear();
        measurement_id_of_noise_dependency_.clear();

        noise_dwell_time_us_preset_ = 5;

        fixed_noise_bandwidth_ = -1;

        gt_timer_.set_timing_in_destruction(false);
        performTiming_ = false;
    }

    NoiseAdjustGadget::~NoiseAdjustGadget()
    {

    }

    int NoiseAdjustGadget::process_config(ACE_Message_Block* mb)
    {
        boost::shared_ptr<std::string> str = this->get_string_value("workingdirectory");
        if ( !str->empty() )
        {
            noise_dependency_folder_ = *str;
        }
        else
        {
            noise_dependency_folder_ = ACE_OS::getenv("GADGETRON_HOME");
            noise_dependency_folder_.append("/dependencies/noise");
        }
        GADGET_MSG("Folder to store noise dependencies is " << noise_dependency_folder_);

        str = this->get_string_value("noise_dependency_prefix");

        if ( !str->empty() )
        {
            noise_dependency_prefix_ = *str;
        }

        performTiming_ = this->get_bool_value("performTiming");

        noise_dwell_time_us_preset_ = (float)this->get_double_value("noise_dwell_time_us_preset");
        if ( noise_dwell_time_us_preset_ == 0 ) noise_dwell_time_us_preset_ = 5;

        boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

        if ( cfg->acquisitionSystemInformation().present() )
        {
            receiver_noise_bandwidth_ = (float)(cfg->acquisitionSystemInformation().get().relativeReceiverNoiseBandwidth().present() ?
                cfg->acquisitionSystemInformation().get().relativeReceiverNoiseBandwidth().get() : 0.793f);

            GADGET_MSG("receiver_noise_bandwidth_ is " << receiver_noise_bandwidth_);
        }

        // find the patient ID
        if ( cfg->subjectInformation().present() )
        {
            if ( cfg->subjectInformation().get().patientID().present() )
            {
                patient_id_ = cfg->subjectInformation().get().patientID().get();
                GADGET_MSG("Patient ID is " << patient_id_);

                size_t len = patient_id_.length();
                for ( size_t n=0; n<len; n++ )
                {
                    if ( patient_id_[n] == '-' ) patient_id_[n] = '_';
                    if ( patient_id_[n] == ':' ) patient_id_[n] = '_';
                }
            }
        }

        // find the study ID
        if ( cfg->studyInformation().present() )
        {
            if ( cfg->studyInformation().get().studyID().present() )
            {
                study_id_ = cfg->studyInformation().get().studyID().get();
                GADGET_MSG("Study ID is " << study_id_);

                size_t len = study_id_.length();
                for ( size_t n=0; n<len; n++ )
                {
                    if ( study_id_[n] == '-' ) study_id_[n] = '_';
                    if ( study_id_[n] == ':' ) study_id_[n] = '_';
                }
            }
        }

        // find the measurementID of this scan
        if ( cfg->measurementInformation().present() )
        {
            if ( cfg->measurementInformation().get().measurementID().present() )
            {
                measurement_id_ = cfg->measurementInformation().get().measurementID().get();
                GADGET_MSG("Measurement ID is " << measurement_id_);
            }

            // find the noise depencies if any
            if ( cfg->measurementInformation().get().measurementDependency().size() > 0 )
            {
                measurement_id_of_noise_dependency_.clear();

                ISMRMRD::measurementInformationType::measurementDependency_sequence& measDependency = cfg->measurementInformation().get().measurementDependency();
                ISMRMRD::measurementInformationType::measurementDependency_sequence::const_iterator iter = measDependency.begin();
                for ( ; iter!=measDependency.end(); iter++ )
                {
                    std::string dependencyType = iter->dependencyType();
                    std::string dependencyID = iter->measurementID();

                    GADGET_MSG("Found dependency measurement : " << dependencyType << " with ID " << dependencyID);

                    if ( dependencyType=="Noise" || dependencyType=="noise" )
                    {
                        measurement_id_of_noise_dependency_ = dependencyID;
                    }
                }

                if ( !measurement_id_of_noise_dependency_.empty() )
                {
                    GADGET_MSG("Measurement ID of noise dependency is " << measurement_id_of_noise_dependency_);

                    full_name_stored_noise_dependency_ = this->generateFullNameWhenLoadNoiseDependency(measurement_id_of_noise_dependency_);
                    GADGET_MSG("Stored noise dependency is " << full_name_stored_noise_dependency_);

                    // try to load the precomputed noise prewhitener
                    if ( !this->loadNoisePrewhitener(noise_dwell_time_us_, noise_covariance_matrixf_) )
                    {
                        GADGET_MSG("Stored noise dependency is NOT found : " << full_name_stored_noise_dependency_);
                        use_stored_noise_prewhitener_ = false;
                        noise_dwell_time_us_ = -1;
                        noise_covariance_matrixf_.clear();
                    }
                    else
                    {
                        GADGET_MSG("Stored noise dependency is found : " << full_name_stored_noise_dependency_);
                        GADGET_MSG("Stored noise dwell time in us is " << noise_dwell_time_us_);
                        GADGET_MSG("Stored noise channel number is " << noise_covariance_matrixf_.get_size(0));
                        use_stored_noise_prewhitener_ = true;
                    }
                }
            }
        }

        // find out the system version if possible
        if ( cfg->userParameters().present() )
        {
            if ( cfg->userParameters().get().userParameterDouble().size() > 0 )
            {
                ISMRMRD::userParameters::userParameterDouble_sequence& userPara = cfg->userParameters().get().userParameterDouble();
                ISMRMRD::userParameters::userParameterDouble_sequence::const_iterator iter = userPara.begin();

                for ( ; iter!=userPara.end(); iter++ )
                {
                    std::string usrParaName = iter->name();
                    double usrParaValue = iter->value();

                    GADGET_MSG("Found user parameter double : " << usrParaName << " with value : " << usrParaValue);

                    if ( usrParaName == "NoiseBandWidth" )
                    {
                        if ( usrParaValue > 0 )
                        {
                            fixed_noise_bandwidth_ = usrParaValue;
                        }
                    }
                }
            }
        }

        return GADGET_OK;
    }

    std::string NoiseAdjustGadget::generateFullNameWhenLoadNoiseDependency(const std::string& measurement_id_of_noise)
    {
        // find the scan prefix
        std::string measurementStr = measurement_id_;
        size_t ind  = measurement_id_.find_last_of ("_");
        if ( ind != std::string::npos )
        {
            measurementStr = measurement_id_.substr(0, ind);
            measurementStr.append("_");
            measurementStr.append(measurement_id_of_noise);
        }

        std::string full_name_loaded_noise_dependency;

        full_name_loaded_noise_dependency = noise_dependency_folder_;
        full_name_loaded_noise_dependency.append("/");
        full_name_loaded_noise_dependency.append(noise_dependency_prefix_);
        full_name_loaded_noise_dependency.append("_");
        full_name_loaded_noise_dependency.append(measurementStr);

        return full_name_loaded_noise_dependency;
    }

    std::string NoiseAdjustGadget::generateFullNameWhenStoreNoiseDependency(const std::string& measurement_id)
    {
        std::string full_name_stored_noise_dependency;

        full_name_stored_noise_dependency = noise_dependency_folder_;
        full_name_stored_noise_dependency.append("/");
        full_name_stored_noise_dependency.append(noise_dependency_prefix_);
        full_name_stored_noise_dependency.append("_");
        full_name_stored_noise_dependency.append(measurement_id);

        return full_name_stored_noise_dependency;
    }

    bool NoiseAdjustGadget::loadNoisePrewhitener(float& noise_dwell_time_us, hoNDArray< ValueType >& noise_covariance_matrixf)
    {
        std::ifstream infile;
        infile.open (full_name_stored_noise_dependency_.c_str(), std::ios::in|std::ios::binary);

        if (infile.good() )
        {
            infile.read( reinterpret_cast<char*>(&noise_dwell_time_us), sizeof(float));

            size_t len;
            infile.read( reinterpret_cast<char*>(&len), sizeof(size_t));

            char* buf = new char[len];
            if ( buf == NULL ) return false;

            infile.read(buf, len);

            if ( !noise_covariance_matrixf.deserialize(buf, len) )
            {
                delete [] buf;
                return false;
            }

            delete [] buf;
            infile.close();
        }
        else
        {
            GADGET_ERROR_MSG("Noise prewhitener file is not good for writing");
            return false;
        }

        return true;
    }

    bool NoiseAdjustGadget::saveNoisePrewhitener(const std::string& full_name_stored_noise_dependency, float& noise_dwell_time_us, hoNDArray< ValueType >& noise_covariance_matrixf)
    {
        char* buf = NULL;
        size_t len(0);
        if ( !noise_covariance_matrixf.serialize(buf, len) )
        {
            GADGET_ERROR_MSG("Noise prewhitener serialization failed ...");
            return false;
        }

        std::ofstream outfile;
        outfile.open (full_name_stored_noise_dependency.c_str(), std::ios::out|std::ios::binary);

        if (outfile.good())
        {
            GADGET_MSG("write out the noise dependency file : " << full_name_stored_noise_dependency);
            outfile.write( reinterpret_cast<char*>(&noise_dwell_time_us), sizeof(float));
            outfile.write( reinterpret_cast<char*>(&len), sizeof(size_t));
            outfile.write(buf, len);
            outfile.close();
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

    void NoiseAdjustGadget::computeNoisePrewhitener(bool savePrewhitener)
    {
        GADGET_START_TIMING_CONDITION(gt_timer_, "compute noise prewhitener ... ", performTiming_);

        if ( noise_dwell_time_us_ > 0 )
        {
            if (number_of_noise_samples_ > 1)
            {
                GADGET_MSG("Noise dwell time: " << noise_dwell_time_us_);
                GADGET_MSG("receiver_noise_bandwidth: " << receiver_noise_bandwidth_);

                if (!noise_decorrelation_calculated_)
                {
                    GADGET_MSG("Calculating noise decorrelation");

                    // Armadillo can best do its template magic when we concatenate all the operations...
                    // 1. scale for number of samples
                    // 2. Cholesky decomposition
                    // 3. Invert lower triangular

                    arma::cx_fmat noise_covf = as_arma_matrix(&noise_covariance_matrixf_);

                    {
                        noise_covf = arma::inv(arma::trimatu(arma::chol(noise_covf/( (float)number_of_noise_samples_-1))));
                    }

                    // save the noise prewhitener
                    if ( savePrewhitener )
                    {
                        std::string fullNameOfStoredNoiseDependency;
                        fullNameOfStoredNoiseDependency = this->generateFullNameWhenStoreNoiseDependency(measurement_id_);
                        this->saveNoisePrewhitener(fullNameOfStoredNoiseDependency, noise_dwell_time_us_, noise_covariance_matrixf_);
                    }

                    noise_decorrelation_calculated_ = true;
                }
            }
        }

        GADGET_STOP_TIMING_CONDITION(gt_timer_, performTiming_);
    }

    int NoiseAdjustGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
    {
        // GADGET_START_TIMING_CONDITION(gt_timer_, "in noise process ... ", performTiming_);

        bool is_scc_correction = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA).isSet(m1->getObjectPtr()->flags);
        bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);

        unsigned int channels = m1->getObjectPtr()->active_channels;
        unsigned int samples = m1->getObjectPtr()->number_of_samples;

        if ( measurement_id_.empty() )
        {
            unsigned int muid = m1->getObjectPtr()->measurement_uid;
            std::ostringstream ostr;
            ostr << muid;
            measurement_id_ = ostr.str();
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
                    arma::cx_fmat noise_covf = as_arma_matrix(&noise_covariance_matrixf_);
                    noise_covf *= noise_bw_scale_factor_;
                    noise_decorrelation_calculated_ = true;
                }

                if (noise_decorrelation_calculated_)
                {
                    // GADGET_START_TIMING_CONDITION(gt_timer_, "apply noise prewhitener ... ", performTiming_);

                    //#ifdef USE_MKL
                    //    GeneralMatrixProduct_gemm(*m2->getObjectPtr(), *m2->getObjectPtr(), false, noise_covariance_matrixf_, false);
                    //#else
                        arma::cx_fmat noise_covf = as_arma_matrix(&noise_covariance_matrixf_);
                        arma::cx_fmat am2 = as_arma_matrix(m2->getObjectPtr());
                        am2 = am2*arma::trimatu(noise_covf);
                    //#endif // USE_MKL

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

                #ifdef USE_MKL
                    GADGET_CHECK_RETURN(GeneralMatrixProduct_gemm(noise_covariance_matrixf_once_, *m2->getObjectPtr(), true, *m2->getObjectPtr(), false), GADGET_FAIL);
                    GADGET_CHECK_RETURN(Gadgetron::add(noise_covariance_matrixf_once_, noise_covariance_matrixf_, noise_covariance_matrixf_), GADGET_FAIL);
                #else
                    for (unsigned int s = 0; s < samples; s++)
                    {
                        for (unsigned int i = 0; i < channels; i++)
                        {
                            for (unsigned int j = 0; j < channels; j++)
                            {
                                cc_ptr[i*channels + j] += (data_ptr[i * samples + s] * conj(data_ptr[j * samples + s]));
                            }
                        }
                    }
                #endif // USE_MKL

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
                                this->computeNoisePrewhitener(true);
                            }
                            else
                            {
                                this->computeNoisePrewhitener(false);
                            }

                            // apply the scaling
                            arma::cx_fmat noise_covf = as_arma_matrix(&noise_covariance_matrixf_);
                            noise_covf *= noise_bw_scale_factor_;
                        }
                        else
                        {
                            if ( m2->getObjectPtr()->get_size(1) == noise_covariance_matrixf_.get_size(0) )
                            {
                                //#ifdef USE_MKL
                                //    GeneralMatrixProduct_gemm(*m2->getObjectPtr(), *m2->getObjectPtr(), false, noise_covariance_matrixf_, false);
                                //#else
                                    arma::cx_fmat noise_covf = as_arma_matrix(&noise_covariance_matrixf_);
                                    arma::cx_fmat am2 = as_arma_matrix(m2->getObjectPtr());
                                    am2 = am2*arma::trimatu(noise_covf);
                                //#endif // USE_MKL
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

    int NoiseAdjustGadget::close(unsigned long flags)
    {
        if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

        if ( !computed_in_close_ )
        {
            computed_in_close_ = true;
            if ( !this->use_stored_noise_prewhitener_ )
            {
                if ( noise_dwell_time_us_ < 0 ) noise_dwell_time_us_ = noise_dwell_time_us_preset_; // this scan is a noise measurement
                this->computeNoisePrewhitener(true);
            }
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(NoiseAdjustGadget)

} // namespace Gadgetron
