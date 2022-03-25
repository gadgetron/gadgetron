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
                    
        std::chrono::milliseconds sleep_time;
        NODE_PROPERTY(sleep_time_, int, "sleep_time", 0);

    };
}
#endif //ACCUMULATORGADGET_H