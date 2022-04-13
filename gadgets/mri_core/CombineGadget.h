/**
    \brief  Combines coils and outputs magnitude images
    \author Original: Souheil Inati
    \author PureGadget Conversion: Andrew Dupuis
    \test   EPI_2d.cfg
*/

#pragma once

#include "PureGadget.h"
#include "Types.h"

namespace Gadgetron{
    class CombineGadget : public Core::PureGadget<Core::AnyImage, Core::AnyImage> {
    public:
      using Core::PureGadget<Core::AnyImage,Core::AnyImage>::PureGadget;
        Core::AnyImage process_function(Core::AnyImage image) const override;
    };
}
