#pragma once

#include "Node.h"
#include "gadgetron_mricore_export.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE ImageFinishGadget : public Core::ChannelGadget {
    public:
        ImageFinishGadget(
                const Core::Context &context,
                const Core::GadgetProperties &properties
        ) : ChannelGadget(properties) {};

    protected:
        void process(Core::InputChannel& in,
                    Core::OutputChannel& out) override;

    };
}

