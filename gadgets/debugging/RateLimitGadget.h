#pragma once

#include "Gadget.h"
#include "hoNDArray.h"


namespace Gadgetron{

    class RateLimitGadget : public Core::ChannelGadget<mrd::StreamItem>
    {
    public:
      using Core::ChannelGadget<mrd::StreamItem>::ChannelGadget;

      RateLimitGadget(const Core::Context& context, const Core::GadgetProperties& props);

      void process(Core::InputChannel<mrd::StreamItem>& input, Core::OutputChannel& output) override;

    protected:
      NODE_PROPERTY(sleep_time_int, int, "sleep_time", 0);

      std::chrono::milliseconds sleep_time_;
    };
}
