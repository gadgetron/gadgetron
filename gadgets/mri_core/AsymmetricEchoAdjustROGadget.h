/**
    \brief  Re-aligns the readout data with the center of the echo at the center of the incoming array (if not a noise scan and partial fourier along readout is detected)
    \author Original: Hui Xue
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Tested by: gpu_grappa_simple.cfg, cpu_grappa_simple.cfg, generic_cartesian_cine_denoise.cfg, and others 
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"

namespace Gadgetron{
  class AsymmetricEchoAdjustROGadget : public Core::ChannelGadget<Core::Acquisition> 
    {
      public:
        using Core::ChannelGadget<Core::Acquisition>::ChannelGadget;
        AsymmetricEchoAdjustROGadget(const Core::Context& context, const Core::GadgetProperties& props);
        ~AsymmetricEchoAdjustROGadget() override = default;
        void process(Core::InputChannel<Core::Acquisition>& input, Core::OutputChannel& output) override;
      protected:
        std::vector<unsigned int> maxRO_;
    };
}
