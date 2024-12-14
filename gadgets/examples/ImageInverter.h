#pragma once

#include "PureGadget.h"

namespace Gadgetron::Examples {
    class ImageInverter : public Core::PureGadget<mrd::AnyImage, mrd::AnyImage> {
    public:
      using Core::PureGadget<mrd::AnyImage, mrd::AnyImage>::PureGadget;
        mrd::AnyImage process_function(mrd::AnyImage image) const override;
    };
}
