#pragma once

#include "Node.h"
#include "gadgetron_mricore_export.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE ImageFinishGadget : public Core::GadgetNode {
    public:
        ImageFinishGadget(const Core::Context &context, const Core::GadgetProperties &properties) : GadgetNode(
                properties) {};

    protected:
        void process(std::shared_ptr<Core::InputChannel<Core::Message>> in,
                     std::shared_ptr<Core::OutputChannel> out) override;

    };
}

