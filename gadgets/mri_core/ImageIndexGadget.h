#pragma once

#include "Node.h"
#include "Types.h"

namespace Gadgetron {

    class ImageIndexGadget : public Core::ChannelGadget<Core::AnyImage> {
      public:
        ImageIndexGadget(const Core::Context &, const Core::GadgetProperties &);
        void process(Core::InputChannel<Core::AnyImage> &, Core::OutputChannel &) override;
    };
}
