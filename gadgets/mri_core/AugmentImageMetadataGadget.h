/** \file   AugmentImageMetadataGadget.h
    \brief  This Gadget augments the image metadata headers by adding some fields that can be required by clients
    \author Michael Hansen
*/

#pragma once
#include "hoNDArray.h"
#include "ismrmrd/meta.h"

#include <ismrmrd/ismrmrd.h>

#include <Types.h>
#include "PureGadget.h"
namespace Gadgetron
{
class AugmentImageMetadataGadget: public Core::PureGadget<Core::Image<std::complex<float>>,Core::Image<std::complex<float>>>
    {
    public:
        AugmentImageMetadataGadget(const Core::Context& context, const Core::GadgetProperties& props);
        Core::Image<std::complex<float>> process_function(Core::Image<std::complex<float>> args) const override;
};
}

