#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"

namespace Gadgetron{
  class AcquisitionPassthroughGadget :
    public Core::ChannelGadget<Core::Acquisition> 
    {
      public:
        AcquisitionPassthroughGadget(const Core::Context& context, const Core::GadgetProperties& props);
        ~AcquisitionPassthroughGadget() override = default;
        void process(Core::InputChannel<Core::Acquisition>& input, Core::OutputChannel& output) override;
    };
}