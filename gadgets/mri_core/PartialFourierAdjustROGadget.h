/**
    \brief  
    \author Original: Hui Xue
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
  class PartialFourierAdjustROGadget : public Core::ChannelGadget<Core::Acquisition> 
    {
      public:
        using Core::ChannelGadget<Core::Acquisition>::ChannelGadget;
        PartialFourierAdjustROGadget(const Core::Context& context, const Core::GadgetProperties& props);
        ~PartialFourierAdjustROGadget() override = default;
        void process(Core::InputChannel<Core::Acquisition>& input, Core::OutputChannel& output) override;
      protected:
        std::vector<unsigned int> maxRO_;
    };
}
