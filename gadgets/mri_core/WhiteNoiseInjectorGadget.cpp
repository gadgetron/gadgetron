#include "WhiteNoiseInjectorGadget.h"
#include "PureGadget.h"
#include "Types.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_math.h"
#include "ismrmrd/xml.h"
#include <array>

namespace Gadgetron {

template <typename T> RandNormGenerator<T>::RandNormGenerator() {
    rng_.seed();
    this->setPara(0, 1);
}

template <typename T> RandNormGenerator<T>::RandNormGenerator(long long s, T mean, T sigma) {
    this->seed(s);
    this->setPara(mean, sigma);
}

template <typename T> RandNormGenerator<T>::~RandNormGenerator() {}

template <typename T> void RandNormGenerator<T>::seed(unsigned long s) { rng_.seed(s); }

template <typename T> void RandNormGenerator<T>::setPara(T mean, T sigma) {
    typename std::normal_distribution<T>::param_type para(mean, sigma);
    dist_norm_.param(para);
}

template <typename T> inline void RandNormGenerator<T>::gen(hoNDArray<T>& randNum) {
    try {
        size_t N = randNum.get_number_of_elements();
        size_t n;
        for (n = 0; n < N; n++) {
            randNum(n) = dist_norm_(rng_);
        }
    } catch (...) {
        GADGET_THROW("Errors in RandNormGenerator<T>::gen(hoNDArray<T>& randNum) ... ");
    }
}

template <typename T> inline void RandNormGenerator<T>::gen(hoNDArray<std::complex<T>>& randNum) {
    try {
        size_t N = randNum.get_number_of_elements();
        size_t n;

        T real, imag;
        for (n = 0; n < N; n++) {
            real = dist_norm_(rng_);
            imag = dist_norm_(rng_);

            randNum(n) = std::complex<T>(real, imag);
        }
    } catch (...) {
        GADGET_THROW("Errors in RandNormGenerator<T>::gen(hoNDArray< std::complex<T> >& randNum) ... ");
    }
}

WhiteNoiseInjectorGadget::~WhiteNoiseInjectorGadget() { delete randn_; }

WhiteNoiseInjectorGadget::WhiteNoiseInjectorGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : Core::ChannelGadget<Core::Acquisition>(context, props) {

    randn_ = new RandNormGenerator<double>();

    acceFactorE1_ = 1;
    acceFactorE2_ = 1;

    is_interleaved_ = false;
    is_embeded_ = false;
    is_seperate_ = false;
    is_external_ = false;
    is_other_ = false;
    is_no_acceleration_ = false;

    GDEBUG_STREAM("noise mean is " << noise_mean);
    GDEBUG_STREAM("noise std is " << noise_std);
    GDEBUG_STREAM("add_noise_ref is " << add_noise_ref);

    randn_->setPara(noise_mean, noise_std);

    // get the current time and generate a seed
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    long long seed =
        (long long)(1e10 * (timeinfo->tm_year + 1900) + 1e8 * (timeinfo->tm_mon + 1) + 1e6 * timeinfo->tm_mday +
                    1e4 * timeinfo->tm_hour + 1e2 * timeinfo->tm_min + timeinfo->tm_sec + std::rand());

    std::array<unsigned int, 10> sequence;
    sequence[0] = (unsigned int)(1e10 * (timeinfo->tm_year + 1900));
    sequence[1] = (unsigned int)(1e8 * (timeinfo->tm_mon + 1));
    sequence[2] = (unsigned int)(1e6 * timeinfo->tm_mday);
    sequence[3] = (unsigned int)(1e4 * timeinfo->tm_hour);
    sequence[4] = (unsigned int)(1e2 * timeinfo->tm_min);
    sequence[5] = (unsigned int)(timeinfo->tm_sec);

    std::srand((unsigned int)seed);
    sequence[6] = (unsigned int)(std::rand());
    sequence[7] = (unsigned int)(std::rand());
    sequence[8] = (unsigned int)(std::rand());
    sequence[9] = (unsigned int)(std::rand());

    std::seed_seq seedSeq(sequence.begin(), sequence.end());
    randn_->getRandomer().seed(seedSeq);

    randn_->seed((unsigned long)seed);

    // ---------------------------------------------------------------------------------------------------------
    ISMRMRD::IsmrmrdHeader h = context.header;

    if (h.encoding.size() != 1) {
        GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
        GDEBUG("This simple WhiteNoiseInjectorGadget only supports one encoding space\n");
        // TODO: How to throw Gadget failures?
    }
    if (!h.encoding[0].parallelImaging) {
        GDEBUG("Parallel Imaging section not found in header");
        // TODO: How to throw Gadget failures?
    }

    ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;

    acceFactorE1_ = (double)(p_imaging.accelerationFactor.kspace_encoding_step_1);
    acceFactorE2_ = (double)(p_imaging.accelerationFactor.kspace_encoding_step_2);

    GDEBUG_STREAM("acceFactorE1_ is " << acceFactorE1_);
    GDEBUG_STREAM("acceFactorE2_ is " << acceFactorE2_);

    if (!p_imaging.calibrationMode.is_present()) {
        GDEBUG("Parallel Imaging calibrationMode not found in header");
        // TODO: How to throw Gadget failures?
    }

    std::string calib = *p_imaging.calibrationMode;
    if (calib.compare("interleaved") == 0) {
        is_interleaved_ = true;
        GDEBUG_STREAM("Calibration mode is interleaved");
    } else if (calib.compare("embedded") == 0) {
        is_embeded_ = true;
        GDEBUG_STREAM("Calibration mode is embedded");
    } else if (calib.compare("separate") == 0) {
        is_seperate_ = true;
        GDEBUG_STREAM("Calibration mode is separate");
    } else if (calib.compare("external") == 0) {
        is_external_ = true;
        GDEBUG_STREAM("Calibration mode is external");
    } else if ((calib.compare("other") == 0)) {
        is_other_ = true;
        GDEBUG_STREAM("Calibration mode is other");
    } else {
        GDEBUG("Failed to process parallel imaging calibration mode");
        // TODO: How to throw Gadget failures?
    }
}

void WhiteNoiseInjectorGadget::process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) {

    for (auto acquisition : in) {
        auto header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);
        auto input_data = std::get<hoNDArray<std::complex<float>>>(acquisition);

        hoNDArray<float> output_data;

        bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(header.flags);
        bool is_scc_correction =
            ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA).isSet(header.flags);
        bool is_ref = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION).isSet(header.flags);
        bool is_ref_kspace =
            ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING).isSet(header.flags);

        size_t channels = header.active_channels;
        size_t samples = header.number_of_samples;

        if (!is_noise && !is_scc_correction) {
            bool add_noise = true;
            if (is_ref && !is_ref_kspace && (is_seperate_ || is_external_)) {
                add_noise = add_noise_ref;
                if (!add_noise) {
                    GDEBUG_STREAM("WhiteNoiseInjectorGadget, noise is not added to the ref acquisitions ... ");
                }
            }
            if (add_noise) {
                if (!noise_.dimensions_equal(&input_data)) {
                    auto dims = input_data.dimensions();
                    noise_.create(dims);
                    noise_fl_.create(dims);
                }

                try {
                    randn_->gen(noise_);
                } catch (...) {
                    GERROR_STREAM("WhiteNoiseInjectorGadget, randn_->gen(noise_) failed ... ");
                    // TODO: How to throw Gadget failures?
                }

                if (!noise_fl_.copyFrom(noise_)) {
                    GERROR_STREAM("WhiteNoiseInjectorGadget, noise_fl_.copyFrom(noise_) failed ... ");
                    // TODO: How to throw Gadget failures?
                }

                try {
                    Gadgetron::add(input_data, noise_fl_, input_data);
                } catch (...) {
                    GERROR_STREAM("WhiteNoiseInjectorGadget, Gadgetron::add(*m2->getObjectPtr(), noise_, "
                                  "*m2->getObjectPtr()) failed ... ");
                    // TODO: How to throw Gadget failures?
                }
            }
        }
        out.push(acquisition);
    }
}

GADGET_FACTORY_DECLARE(WhiteNoiseInjectorGadget)
} // namespace Gadgetron
