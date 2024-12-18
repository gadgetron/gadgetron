#include "RateLimitGadget.h"

#include <thread>
#include <chrono>

namespace Gadgetron {

RateLimitGadget::RateLimitGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : Core::ChannelGadget<mrd::StreamItem>(context, props)
{
  this->sleep_time_ = std::chrono::milliseconds(this->sleep_time_int);
}

void RateLimitGadget::process(Core::InputChannel<mrd::StreamItem>& in, Core::OutputChannel& out)
{
  for (auto item : in) {
    std::this_thread::sleep_for(this->sleep_time_);
    out.push(item);
  }
}

GADGETRON_GADGET_EXPORT(RateLimitGadget)
} // namespace Gadgetron
