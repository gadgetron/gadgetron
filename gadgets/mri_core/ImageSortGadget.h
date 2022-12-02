/**
    \brief  Sorts all pipeline images by the selected sorting dimension flag
    \test   Tested by: distributed_simple_gre.cfg, distributed_buffer_simple_gre.cfg
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"

namespace Gadgetron{

  struct ImageEntry
  {
    Core::AnyImage image_;
    int index_;
  };

  class ImageSortGadget : public Core::ChannelGadget<Core::AnyImage> {
    public:
      GADGET_DECLARE(ImageSortGadget);
      ImageSortGadget(const Core::Context &, const Core::GadgetProperties &);
      void process(Core::InputChannel<Core::AnyImage> &, Core::OutputChannel &) override;
    protected:
      // { "average", "slice", "contrast", "phase", "repetition", "set" }
      NODE_PROPERTY(sorting_dimension, std::string, "Dimension that data will be sorted by", "slice"); 
      std::vector<ImageEntry> images_;
  };
}