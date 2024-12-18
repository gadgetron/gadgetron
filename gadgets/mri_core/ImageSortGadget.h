/**
    \brief  Sorts all pipeline images by the selected sorting dimension flag
    \test   Tested by: distributed_simple_gre.cfg, distributed_buffer_simple_gre.cfg
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "Node.h"

namespace Gadgetron{

  struct ImageEntry
  {
    mrd::AnyImage image_;
    int index_;
  };

  class ImageSortGadget : public Core::ChannelGadget<mrd::AnyImage> {
    public:
      ImageSortGadget(const Core::Context &, const Core::GadgetProperties &);
      void process(Core::InputChannel<mrd::AnyImage> &, Core::OutputChannel &) override;
    protected:
      // { "average", "slice", "contrast", "phase", "repetition", "set" }
      NODE_PROPERTY(sorting_dimension, std::string, "Dimension that data will be sorted by", "slice");
      std::vector<ImageEntry> images_;
  };
}