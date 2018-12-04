#pragma once

#include <Node.h>
#include <mri_core_data.h>

namespace Gadgetron::Core::Gadgets {

    class OversamplingGadget : public TypedGadgetNode<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float>> {



    };
}
