/**
    \brief  
    \author Original: Thomas Sangild Sorensen
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Untested
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"

namespace Gadgetron{
  class NoiseAdjustGadget_unoptimized : public Core::ChannelGadget<Core::Acquisition> 
    {
      public:
        using Core::ChannelGadget<Core::Acquisition>::ChannelGadget;
        NoiseAdjustGadget_unoptimized(const Core::Context& context, const Core::GadgetProperties& props);
        ~NoiseAdjustGadget_unoptimized() override = default;
        void process(Core::InputChannel<Core::Acquisition>& input, Core::OutputChannel& output) override;
      protected:
        bool noise_decorrelation_calculated_;
        hoNDArray<std::complex<double>> noise_covariance_matrix_;
        unsigned long int number_of_noise_samples_;
        float noise_dwell_time_us_;
        float acquisition_dwell_time_us_;
        float noise_bw_scale_factor_;
        float receiver_noise_bandwidth_;
        bool is_configured_;
    };
}
