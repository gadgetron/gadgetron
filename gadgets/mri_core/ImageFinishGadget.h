#pragma once

#include "Node.h"
#include "gadgetron_mricore_export.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE ImageFinishGadget : public Core::GadgetNode {
    public:
        ImageFinishGadget(
                const Core::Context &context,
                const Core::GadgetProperties &properties
        ) : GadgetNode(properties) {};

    protected:
        void process(Core::InputChannel& in,
                    Core::OutputChannel& out) override;

    };
}

