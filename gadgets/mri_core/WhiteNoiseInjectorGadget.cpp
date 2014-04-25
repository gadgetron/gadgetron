
#include "WhiteNoiseInjectorGadget.h"
#include "GadgetIsmrmrdReadWrite.h"

#include "gtPlusUtil.h"

#include <array>

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
    noise_mean_ = this->get_double_value("noise_mean");
    noise_std_ = this->get_double_value("noise_std");
    add_noise_ref_ = this->get_bool_value("add_noise_ref");

    GADGET_MSG("noise mean is " << noise_mean_);
    GADGET_MSG("noise std is " << noise_std_);
    GADGET_MSG("add_noise_ref is " << add_noise_ref_);

    randn_->setPara(noise_mean_, noise_std_);

    // get the current time and generate a seed
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    long long seed = 1e10*(timeinfo->tm_year+1900) + 1e8*(timeinfo->tm_mon+1) + 1e6*timeinfo->tm_mday + 1e4*timeinfo->tm_hour + 1e2*timeinfo->tm_min + timeinfo->tm_sec + std::rand();

    std::array<unsigned int, 10> sequence;
    sequence[0] = 1e10*(timeinfo->tm_year+1900);
    sequence[1] = 1e8*(timeinfo->tm_mon+1);
    sequence[2] = 1e6*timeinfo->tm_mday;
    sequence[3] = 1e4*timeinfo->tm_hour;
    sequence[4] = 1e2*timeinfo->tm_min;
    sequence[5] = timeinfo->tm_sec;

    std::srand(seed);
    sequence[6] = std::rand();
    sequence[7] = std::rand();
    sequence[8] = std::rand();
    sequence[9] = std::rand();

    std::seed_seq seedSeq(sequence.begin(), sequence.end());
    randn_->getRandomer().seed(seedSeq);

    randn_->seed(seed);

// ---------------------------------------------------------------------------------------------------------
    // pass the xml file
    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    // seq object
    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
    if (e_seq.size() != 1)
    {
        GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
        GADGET_DEBUG1("This simple WhiteNoiseInjectorGadget only supports one encoding space\n");
        return GADGET_FAIL;
    }

    // find out the PAT mode
    ISMRMRD::ismrmrdHeader::parallelImaging_optional p_imaging_type = cfg->parallelImaging();
    ISMRMRD::parallelImagingType p_imaging = *p_imaging_type;

    acceFactorE1_ = (size_t)(p_imaging.accelerationFactor().kspace_encoding_step_1());
    acceFactorE2_ = (size_t)(p_imaging.accelerationFactor().kspace_encoding_step_2());
    GADGET_MSG("acceFactorE1_ is " << acceFactorE1_);
    GADGET_MSG("acceFactorE2_ is " << acceFactorE2_);

    ISMRMRD::calibrationModeType calib = *(p_imaging.calibrationMode());
    if ( calib == ISMRMRD::calibrationModeType::interleaved )
    {
        is_interleaved_ = true;
        GADGET_MSG("Calibration mode is interleaved");
    }

    if ( calib == ISMRMRD::calibrationModeType::embedded )
    {
        is_embeded_ = true;
        GADGET_MSG("Calibration mode is embedded");
    }

    if ( calib == ISMRMRD::calibrationModeType::separate )
    {
        is_seperate_ = true;
        GADGET_MSG("Calibration mode is separate");
    }

    if ( calib == ISMRMRD::calibrationModeType::external )
    {
        is_external_ = true;
        GADGET_MSG("Calibration mode is external");
    }

    if ( calib == ISMRMRD::calibrationModeType::other )
    {
        is_other_ = true;
        GADGET_MSG("Calibration mode is other");

        if ( acceFactorE1_==1 && acceFactorE2_==1 )
        {
            is_no_acceleration_ = true;
            GADGET_MSG("No acceleration is used ... ");
        }
    }

    return GADGET_OK;
}

int WhiteNoiseInjectorGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    bool is_scc_correction = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA).isSet(m1->getObjectPtr()->flags);

    bool is_ref = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_PARALLEL_CALIBRATION).isSet(m1->getObjectPtr()->flags);
    bool is_ref_kspace = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING).isSet(m1->getObjectPtr()->flags);

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
                GADGET_MSG("WhiteNoiseInjectorGadget, noise is not added to the ref acquisitions ... ");
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

            if ( !Gadgetron::add(*m2->getObjectPtr(), noise_fl_, *m2->getObjectPtr()) )
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
