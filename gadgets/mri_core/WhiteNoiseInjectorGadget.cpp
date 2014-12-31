#include "WhiteNoiseInjectorGadget.h"
#include "gtPlusUtil.h"
#include <array>
#include "ismrmrd/xml.h"

namespace Gadgetron
{

WhiteNoiseInjectorGadget::WhiteNoiseInjectorGadget() : noise_mean_(0), noise_std_(1.0f)
{
    add_noise_ref_ = true;
    randn_ = new RandGenType();

    acceFactorE1_ = 1;
    acceFactorE2_ = 1;

    is_interleaved_ = false;
    is_embeded_ = false;
    is_seperate_ = false;
    is_external_ = false;
    is_other_ = false;
    is_no_acceleration_ = false;
}

WhiteNoiseInjectorGadget::~WhiteNoiseInjectorGadget()
{
    delete randn_;
}

int WhiteNoiseInjectorGadget::process_config(ACE_Message_Block* mb)
{
    noise_mean_ = (float)this->get_double_value("noise_mean");
    noise_std_ = (float)this->get_double_value("noise_std");
    add_noise_ref_ = this->get_bool_value("add_noise_ref");

    GADGET_MSG_DEPRECATED("noise mean is " << noise_mean_);
    GADGET_MSG_DEPRECATED("noise std is " << noise_std_);
    GADGET_MSG_DEPRECATED("add_noise_ref is " << add_noise_ref_);

    randn_->setPara(noise_mean_, noise_std_);

    // get the current time and generate a seed
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    long long seed = (long long)(1e10*(timeinfo->tm_year+1900) + 1e8*(timeinfo->tm_mon+1) + 1e6*timeinfo->tm_mday + 1e4*timeinfo->tm_hour + 1e2*timeinfo->tm_min + timeinfo->tm_sec + std::rand());

    std::array<unsigned int, 10> sequence;
    sequence[0] = (unsigned int)(1e10*(timeinfo->tm_year+1900));
    sequence[1] = (unsigned int)(1e8*(timeinfo->tm_mon+1));
    sequence[2] = (unsigned int)(1e6*timeinfo->tm_mday);
    sequence[3] = (unsigned int)(1e4*timeinfo->tm_hour);
    sequence[4] = (unsigned int)(1e2*timeinfo->tm_min);
    sequence[5] = (unsigned int)(timeinfo->tm_sec);

    std::srand( (unsigned int)seed );
    sequence[6] = (unsigned int)(std::rand());
    sequence[7] = (unsigned int)(std::rand());
    sequence[8] = (unsigned int)(std::rand());
    sequence[9] = (unsigned int)(std::rand());

    std::seed_seq seedSeq(sequence.begin(), sequence.end());
    randn_->getRandomer().seed(seedSeq);

    randn_->seed( (unsigned long)seed );

// ---------------------------------------------------------------------------------------------------------
    ISMRMRD::IsmrmrdHeader h;
    try {
      deserialize(mb->rd_ptr(),h);
    } catch (...) {
      GDEBUG("Error parsing ISMRMRD Header");
      throw;
      return GADGET_FAIL;
    }

    if( h.encoding.size() != 1)
    {
      GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
      GDEBUG("This simple WhiteNoiseInjectorGadget only supports one encoding space\n");
      return GADGET_FAIL;
    }
    if (!h.encoding[0].parallelImaging) {
      GDEBUG("Parallel Imaging section not found in header");
      return GADGET_FAIL;
    }

    ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;

    acceFactorE1_ = (double)(p_imaging.accelerationFactor.kspace_encoding_step_1);
    acceFactorE2_ = (double)(p_imaging.accelerationFactor.kspace_encoding_step_2);

    GADGET_MSG_DEPRECATED("acceFactorE1_ is " << acceFactorE1_);
    GADGET_MSG_DEPRECATED("acceFactorE2_ is " << acceFactorE2_);

    if ( !p_imaging.calibrationMode.is_present() )
    {
        GDEBUG("Parallel Imaging calibrationMode not found in header");
        return GADGET_FAIL;
    }

    std::string calib = *p_imaging.calibrationMode;
    if ( calib.compare("interleaved") == 0 )
    {
      is_interleaved_ = true;
      GADGET_MSG_DEPRECATED("Calibration mode is interleaved");
    } else if ( calib.compare("embedded") == 0 ) {
      is_embeded_ = true;
      GADGET_MSG_DEPRECATED("Calibration mode is embedded");
    } else if ( calib.compare("separate") == 0 ) {
      is_seperate_ = true;
      GADGET_MSG_DEPRECATED("Calibration mode is separate");
    } else if ( calib.compare("external") == 0 ) {
      is_external_ = true;
      GADGET_MSG_DEPRECATED("Calibration mode is external");
    } else if ( (calib.compare("other") == 0)) {
      is_other_ = true;
      GADGET_MSG_DEPRECATED("Calibration mode is other");
    } else {
      GDEBUG("Failed to process parallel imaging calibration mode");
      return GADGET_FAIL;
    }
    
    return GADGET_OK;
}

int WhiteNoiseInjectorGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    bool is_scc_correction = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA).isSet(m1->getObjectPtr()->flags);

    bool is_ref = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION).isSet(m1->getObjectPtr()->flags);
    bool is_ref_kspace = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING).isSet(m1->getObjectPtr()->flags);

    size_t channels = m1->getObjectPtr()->active_channels;
    size_t samples = m1->getObjectPtr()->number_of_samples;

    if (!is_noise && !is_scc_correction )
    {
        bool add_noise = true;
        if ( is_ref && !is_ref_kspace && (is_seperate_||is_external_) )
        {
            add_noise = add_noise_ref_;

            if ( !add_noise )
            {
                GADGET_MSG_DEPRECATED("WhiteNoiseInjectorGadget, noise is not added to the ref acquisitions ... ");
            }
        }

        if ( add_noise )
        {
            if ( !noise_.dimensions_equal(m2->getObjectPtr()) )
            {
                noise_.create(m2->getObjectPtr()->get_dimensions());
                noise_fl_.create(m2->getObjectPtr()->get_dimensions());
            }

            if ( !randn_->gen(noise_) )
            {
                GADGET_ERROR_MSG("WhiteNoiseInjectorGadget, randn_->gen(noise_) failed ... ");
                return GADGET_FAIL;
            }

            if ( !noise_fl_.copyFrom(noise_) )
            {
                GADGET_ERROR_MSG("WhiteNoiseInjectorGadget, noise_fl_.copyFrom(noise_) failed ... ");
                return GADGET_FAIL;
            }

            try
            {
                Gadgetron::add(*m2->getObjectPtr(), noise_fl_, *m2->getObjectPtr());
            }
            catch(...)
            {
                GADGET_ERROR_MSG("WhiteNoiseInjectorGadget, Gadgetron::add(*m2->getObjectPtr(), noise_, *m2->getObjectPtr()) failed ... ");
                return GADGET_FAIL;
            }
        }
    }

    if (this->next()->putq(m1) == -1) 
    {
        ACE_ERROR_RETURN( (LM_ERROR,
                ACE_TEXT("%p\n"),
                ACE_TEXT("WhiteNoiseInjectorGadget::process, passing data on to next gadget")),
                -1);
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(WhiteNoiseInjectorGadget)
}
