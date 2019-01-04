#pragma once

#include "mri_core_data.h"
#include "hoNDArray.h"

#include "Gadget.h"
#include "Context.h"

namespace Gadgetron::Core::Gadgets {

    class OversamplingGadget : public TypedGadgetNode<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float>>> {
    public:
        OversamplingGadget(const Context &context, std::unordered_map<std::string, std::string> props);

    private:

    };
}
