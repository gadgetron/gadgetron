/**
    \brief  Imposes an artificial ms-delay as a rate limit at a point in the Gadgetron pipeline
    \author Original: David Christoffer Hansen
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Tested by: simple_gre_ratelimit.cfg
*/

#pragma once

#include "Gadget.h"
#include "Node.h"
#include "gadgetron_mricore_export.h"
#include <chrono>
#include <complex>
#include <ismrmrd/ismrmrd.h>

namespace Gadgetron {

    class EXPORTGADGETSMRICORE RateLimitGadget : public Core::GenericChannelGadget {
    public:
        RateLimitGadget(
                const Core::Context &context,
                const Core::GadgetProperties &properties
        ) : GenericChannelGadget(context,properties) {};

    protected:
        void process(Core::GenericInputChannel& in,
                    Core::OutputChannel& out) override;
                    
        NODE_PROPERTY(sleep_time, int, "Sleep Time", 0);

    };
}