#pragma once

#include <ismrmrd/ismrmrd.h>

#include "hoNDArray.h"

#include "Types.h"
#include "Writer.h"

namespace Gadgetron::Core::Writers {

    class AcquisitionWriter :
            public TypedWriter<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float>>, optional<hoNDArray<float>>> {
    protected:
        void serialize(
                std::ostream &stream,
                const ISMRMRD::AcquisitionHeader &header,
                const hoNDArray<std::complex<float>> &data,
                const optional<hoNDArray<float>> &trajectory
        ) override;
    };
}