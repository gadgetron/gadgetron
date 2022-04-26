#include "WhiteNoiseInjectorGadget.h"
#include "hoNDArray_elemwise.h"
#include <array>
#include "ismrmrd/xml.h"

namespace Gadgetron
{

template <typename T>
RandNormGenerator<T>::RandNormGenerator()
{
    rng_.seed();
    this->setPara(0, 1);
}

template <typename T>
RandNormGenerator<T>::RandNormGenerator(long long s, T mean, T sigma)
{
    this->seed(s);
    this->setPara(mean, sigma);
}

template <typename T>
RandNormGenerator<T>::~RandNormGenerator()
{
}

template <typename T>
void RandNormGenerator<T>::seed(unsigned long s)
{
    rng_.seed(s);
}

template <typename T>
void RandNormGenerator<T>::setPara(T mean, T sigma)
{
    typename std::normal_distribution<T>::param_type para(mean, sigma);
    dist_norm_.param(para);
}

template <typename T>
inline void RandNormGenerator<T>::gen(hoNDArray<T>& randNum)
{
    try
    {
        size_t N = randNum.get_number_of_elements();
        size_t n;
        for (n = 0; n<N; n++)
        {
            randNum(n) = dist_norm_(rng_);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in RandNormGenerator<T>::gen(hoNDArray<T>& randNum) ... ");
    }
}

template <typename T>
inline void RandNormGenerator<T>::gen(hoNDArray< std::complex<T> >& randNum)
{
    try
    {
        size_t N = randNum.get_number_of_elements();
        size_t n;

        T real, imag;
        for (n = 0; n<N; n++)
        {
            real = dist_norm_(rng_);
            imag = dist_norm_(rng_);

            randNum(n) = std::complex<T>(real, imag);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in RandNormGenerator<T>::gen(hoNDArray< std::complex<T> >& randNum) ... ");
    }
}

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
    noise_mean_ = noise_mean.value();
    noise_std_ = noise_std.value();
    add_noise_ref_ = add_noise_ref.value();

    GDEBUG_STREAM("noise mean is " << noise_mean_);
    GDEBUG_STREAM("noise std is " << noise_std_);
    GDEBUG_STREAM("add_noise_ref is " << add_noise_ref_);

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

    GDEBUG_STREAM("acceFactorE1_ is " << acceFactorE1_);
    GDEBUG_STREAM("acceFactorE2_ is " << acceFactorE2_);

    if ( !p_imaging.calibrationMode.is_present() )
    {
        GDEBUG("Parallel Imaging calibrationMode not found in header");
        return GADGET_FAIL;
    }

    std::string calib = *p_imaging.calibrationMode;
    if ( calib.compare("interleaved") == 0 )
    {
      is_interleaved_ = true;
      GDEBUG_STREAM("Calibration mode is interleaved");
    } else if ( calib.compare("embedded") == 0 ) {
      is_embeded_ = true;
      GDEBUG_STREAM("Calibration mode is embedded");
    } else if ( calib.compare("separate") == 0 ) {
      is_seperate_ = true;
      GDEBUG_STREAM("Calibration mode is separate");
    } else if ( calib.compare("external") == 0 ) {
      is_external_ = true;
      GDEBUG_STREAM("Calibration mode is external");
    } else if ( (calib.compare("other") == 0)) {
      is_other_ = true;
      GDEBUG_STREAM("Calibration mode is other");
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
                GDEBUG_STREAM("WhiteNoiseInjectorGadget, noise is not added to the ref acquisitions ... ");
            }
        }

        if ( add_noise )
        {
            if ( !noise_.dimensions_equal(m2->getObjectPtr()) )
            {
                noise_.create(m2->getObjectPtr()->dimensions());
                noise_fl_.create(m2->getObjectPtr()->dimensions());
            }

            try
            {
                randn_->gen(noise_);
            }
            catch(...)
            {
                GERROR_STREAM("WhiteNoiseInjectorGadget, randn_->gen(noise_) failed ... ");
                return GADGET_FAIL;
            }

            if ( !noise_fl_.copyFrom(noise_) )
            {
                GERROR_STREAM("WhiteNoiseInjectorGadget, noise_fl_.copyFrom(noise_) failed ... ");
                return GADGET_FAIL;
            }

            try
            {
                Gadgetron::add(*m2->getObjectPtr(), noise_fl_, *m2->getObjectPtr());
            }
            catch(...)
            {
                GERROR_STREAM("WhiteNoiseInjectorGadget, Gadgetron::add(*m2->getObjectPtr(), noise_, *m2->getObjectPtr()) failed ... ");
                return GADGET_FAIL;
            }
        }
    }

    if (this->next()->putq(m1) == -1) 
    {
      GERROR("WhiteNoiseInjectorGadget::process, passing data on to next gadget");
      return -1;
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(WhiteNoiseInjectorGadget)
}
