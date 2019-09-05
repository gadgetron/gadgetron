#pragma once

#include "PureGadget.h"
#include "Types.h"

namespace Gadgetron::Examples {
    class ImageInverter : public Core::PureGadget<Core::AnyImage, Core::AnyImage> {
    public:
        ImageInverter(const Core::Context &, const Core::GadgetProperties &);
        Core::AnyImage process_function(Core::AnyImage image) const override;
    };
}
