#include "RateLimitGadget.h"
#include "ismrmrd/xml.h"
#include <thread>
#include <chrono>

#include "ImageFinishGadget.h"

namespace Gadgetron {

    void RateLimitGadget::process(Core::GenericInputChannel& in, Core::OutputChannel& out) {
        for (auto message : in) {
            std::string msg = "Rate Limit || Sleep for: " + std::to_string(sleep_time) + "ms\n";
            GDEBUG(msg.c_str());
            std::this_thread::sleep_for(std::chrono::milliseconds(sleep_time));
            out.push_message(std::move(message));
            GDEBUG("Rate Limit || Finished Sleep\n");
        }
    }

    GADGETRON_GADGET_EXPORT(RateLimitGadget); // TODO: difference from GADGET_FACTORY_DECLARE(RateLimitGadget)?
}
