#pragma once

#include <ismrmrd/ismrmrd.h>

#include "hoNDArray.h"

#include "Types.h"
#include "Writer.h"

namespace Gadgetron::Core::Writers {

    class AcquisitionWriter :
            public TypedWriter<ISMRMRD::AcquisitionHeader, optional<hoNDArray<float>>, hoNDArray<std::complex<float>>> {
    protected:
        void serialize(
                std::ostream &stream,
                const ISMRMRD::AcquisitionHeader &header,
                const optional<hoNDArray<float>> &trajectory,
                const hoNDArray<std::complex<float>> &data
        ) override;
    };
}