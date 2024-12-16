#pragma once

#include "Node.h"

namespace Gadgetron {

    /**
     * ImageIndexGadget sets the image_index fields for images.
     *
     * ImageIndexGadget will keep track of individual image series, and number images in each series
     * sequentially.
     */
    class ImageIndexGadget : public Core::ChannelGadget<mrd::AnyImage> {
      public:
        ImageIndexGadget(const Core::Context &, const Core::GadgetProperties &);
        void process(Core::InputChannel<mrd::AnyImage> &, Core::OutputChannel &) override;
    };
}
