#ifndef RATELIMITGADGET_H
#define RATELIMITGADGET_H

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
#endif //ACCUMULATORGADGET_H