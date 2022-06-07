/**
    \brief  Combines coils and outputs magnitude images
    \test   EPI_2d.cfg
*/

#pragma once

#include "PureGadget.h"
#include "Types.h"
#include "hoNDArray_math.h"

namespace Gadgetron{
    class CombineGadget : public Core::PureGadget<Core::AnyImage, Core::AnyImage> {
    public:
      using Core::PureGadget<Core::AnyImage,Core::AnyImage>::PureGadget;
        Core::AnyImage process_function(Core::AnyImage image) const override;
    };
}
