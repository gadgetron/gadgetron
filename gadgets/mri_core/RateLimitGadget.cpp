#include "RateLimitGadget.h"
#include "ismrmrd/xml.h"
#include <thread>
#include <chrono>

#include "ImageFinishGadget.h"

namespace Gadgetron {

    void RateLimitGadget::process(Core::GenericInputChannel& in, Core::OutputChannel& out) {
        for (auto message : in) {
            this->sleep_time = std::chrono::milliseconds(sleep_time_);
            std::this_thread::sleep_for(this->sleep_time);
            out.push_message(std::move(message));
        }
    }

    GADGETRON_GADGET_EXPORT(RateLimitGadget); // TODO: difference from GADGET_FACTORY_DECLARE(RateLimitGadget)?
}
