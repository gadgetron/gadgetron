/**
    \brief  Adds parameterized white noise to an incoming acquisition
    \author Original: Hui Xue
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Untested
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "gadgetron_mricore_export.h"
#include <random>
#include "PureGadget.h"

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


    class WhiteNoiseInjectorGadget : public Core::ChannelGadget<Core::Acquisition> {
    public:
        WhiteNoiseInjectorGadget(const Core::Context& context, const Core::GadgetProperties& props);
        ~WhiteNoiseInjectorGadget() override;
        void process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) override;
    protected:
        NODE_PROPERTY(noise_mean, float, "Noise mean", 0.0);
        NODE_PROPERTY(noise_std, float, "Noise standard deviation", 0.0);
        NODE_PROPERTY(add_noise_ref, bool, "Add noise to reference scans", false);

        /// random noise generator
        RandNormGenerator<double>* randn_;

        /// helper memory to store noise
        hoNDArray<std::complex<double>> noise_;
        hoNDArray<std::complex<float>> noise_fl_;

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

