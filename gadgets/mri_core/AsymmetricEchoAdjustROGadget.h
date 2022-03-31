#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"

/// for incoming readout
/// if not the noise scan and the partial fourier along readout is detected
/// the readout data will be realigned with center of echo at the centre of incoming 1D array
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
