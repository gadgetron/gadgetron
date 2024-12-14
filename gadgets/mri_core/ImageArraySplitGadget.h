/**
    \brief  Splits an ImageArray and outputs separate images
    \test   simple_gre_3d.cfg, distributed_simple_gre.cfg, and others
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "Node.h"
#include "hoNDArray_math.h"

namespace Gadgetron{

  using ImageOrImageArray = std::variant<mrd::AnyImage, mrd::ImageArray>;

  class ImageArraySplitGadget : public Core::ChannelGadget<ImageOrImageArray>
    {
      public:
        using Core::ChannelGadget<ImageOrImageArray>::ChannelGadget;
        void process(Core::InputChannel<ImageOrImageArray>& input, Core::OutputChannel& output) override;
    };
}
