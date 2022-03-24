#pragma once

#ifndef IMAGESORTGADGET_H
#define IMAGESORTGADGET_H

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
      //GadgetPropertyLimitsEnumeration<std::string> limits { "average", "slice", "contrast", "phase", "repetition", "set" };
      NODE_PROPERTY(sorting_dimension, std::string, "Dimension that data will be sorted by", "slice"); // TODO: Add property limits back on? 
      std::vector<ImageEntry> images_;
  };
}

#endif //IMAGESORTGADGET_H
