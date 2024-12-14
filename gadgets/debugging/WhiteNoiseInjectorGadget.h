
#pragma once

#include "Gadget.h"
#include "hoNDArray.h"

#include <random>

namespace Gadgetron
{

template <typename T>
class RandNormGenerator
{
public:

    typedef std::mt19937 RandomGeneratorType;

    RandNormGenerator();
    RandNormGenerator(long long seed, T mean = 0, T sigma = 1);
    ~RandNormGenerator();

    void seed(unsigned long seed);
    void setPara(T mean = 0, T sigma = 1);

    RandomGeneratorType& getRandomer() { return rng_; }
    const RandomGeneratorType& getRandomer() const { return rng_; }

    void gen(hoNDArray<T>& randNum);
    void gen(hoNDArray< std::complex<T> >& randNum);

protected:

    RandomGeneratorType rng_;
    std::normal_distribution<T> dist_norm_;
};

/// add white noise to the kspace data
class WhiteNoiseInjectorGadget : public Gadgetron::Gadget1<mrd::Acquisition>
{
public:
    typedef Gadgetron::RandNormGenerator<double> RandGenType;

    WhiteNoiseInjectorGadget();
    virtual ~WhiteNoiseInjectorGadget();

protected:
    GADGET_PROPERTY(noise_mean, float, "Noise mean", 0.0);
    GADGET_PROPERTY(noise_std, float, "Noise standard deviation", 0.0);
    GADGET_PROPERTY(add_noise_ref, bool, "Add noise to reference scans", false);

    virtual int process_config(const mrd::Header& header);

    virtual int process(Gadgetron::GadgetContainerMessage<mrd::Acquisition>* m1);

    /// whether to add noise to ref acquisition
    bool add_noise_ref_;

    /// noise mean and standard deviation
    float noise_mean_;
    float noise_std_;

    /// random noise generator
    RandGenType* randn_;

    /// helper memory to store noise
    hoNDArray< std::complex<double> > noise_;
    hoNDArray< std::complex<float> > noise_fl_;

    /// calibration mode and rate
    size_t acceFactorE1_;
    size_t acceFactorE2_;

    bool is_interleaved_;
    bool is_embeded_;
    bool is_seperate_;
    bool is_external_;
    bool is_other_;
    bool is_no_acceleration_;
};

}
