/**
    \brief  Resizes images to either a deisred size or by a desired scaling factor
    \author Original: Hui Xue
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Untested
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"

namespace Gadgetron
{
  class ImageResizingGadget : public Core::ChannelGadget<Core::AnyImage> {
    public:
      GADGET_DECLARE(ImageResizingGadget);
      ImageResizingGadget(const Core::Context &, const Core::GadgetProperties &);
      void process(Core::InputChannel<Core::AnyImage> &, Core::OutputChannel &) override;
    protected:
      NODE_PROPERTY(new_RO, size_t, "New image size along RO; if 0, no effect", 0);
      NODE_PROPERTY(new_E1, size_t, "New image size along E1; if 0, no effect", 0);
      NODE_PROPERTY(new_E2, size_t, "New image size along E2; if 0, no effect", 0);

      NODE_PROPERTY(scale_factor_RO, double, "Scale factors; if 0, no effect", 1.0);
      NODE_PROPERTY(scale_factor_E1, double, "Scale factors; if 0, no effect", 1.0);
      NODE_PROPERTY(scale_factor_E2, double, "Scale factors; if 0, no effect", 1.0);

      // BSpline interpolation was used
      NODE_PROPERTY(order_interpolator, size_t, "Order of interpolator; higher order may increase noise level", 5);
  };
}
