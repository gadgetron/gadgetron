/**
    \brief  Combines coils and outputs magnitude images
    \test   EPI_2d.cfg
*/

#pragma once

#include "PureGadget.h"
#include "hoNDArray_math.h"

namespace Gadgetron{
    class CombineGadget : public Core::PureGadget<mrd::AnyImage, mrd::AnyImage> {
        public:
        using Core::PureGadget<mrd::AnyImage, mrd::AnyImage>::PureGadget;
        mrd::AnyImage process_function(mrd::AnyImage image) const override;
    };
}
