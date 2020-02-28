#pragma once

#include <ismrmrd/ismrmrd.h>

#include "hoNDArray.h"

#include "Types.h"
#include "Writer.h"

namespace Gadgetron::Core::Writers {

    class WaveformWriter
            : public Core::TypedWriter<ISMRMRD::WaveformHeader, hoNDArray<uint32_t>> {
    protected:
        void serialize(
                std::ostream &stream,
                const ISMRMRD::WaveformHeader &header,
                const hoNDArray<uint32_t> &array
        ) override;
    };
}