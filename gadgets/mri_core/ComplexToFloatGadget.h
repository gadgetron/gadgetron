/** \file   ComplexToFloatGadget.h
    \brief  This Gadget converts complex float values to float format.
    \author Hui Xue
*/

#pragma once
#include "hoNDArray.h"
#include "ismrmrd/meta.h"

#include <ismrmrd/ismrmrd.h>

#include <Types.h>
#include "PureGadget.h"
namespace Gadgetron
{
class ComplexToFloatGadget: public Core::PureGadget<Core::Image<float>,Core::Image<std::complex<float>>>
    {
    public:
        ComplexToFloatGadget(const Core::Context& context, const Core::GadgetProperties& props);

        Core::Image<float> process_function(Core::Image<std::complex<float>> args) const override;
    private:
        std::map<uint16_t,std::function<hoNDArray<float>(const hoNDArray<std::complex<float>>&)>> converters;
};
}

