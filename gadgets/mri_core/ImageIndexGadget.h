/**
    \brief  Sets the image_index fields for images. Keeps track of individual image series, and number images in each series sequentially
    \author Original: Kristoffer Langeland Knudsen
    \test   Tested by: distributed_simple_gre.cfg
*/

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
