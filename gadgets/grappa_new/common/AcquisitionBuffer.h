#pragma once

#include "Context.h"
#include "Channel.h"
#include "Types.h"

#include "hoNDArray.h"


namespace Gadgetron::Grappa {

    class AcquisitionBuffer {
        using Context = Core::Context;
        using Acquisition = Core::Acquisition;
    public:
        explicit AcquisitionBuffer(Context);

        void add_acquisition(const Acquisition &acquisition);

        hoNDArray<std::complex<float>>
        take(size_t index);

        const hoNDArray<std::complex<float>> &
        view(size_t index);

        void clear(size_t index);

    private:
        const Context context;

        struct {
            std::vector<unsigned long> buffer_dimensions;
            int line_offset;
        } internals;

        std::vector<hoNDArray<std::complex<float>>> buffers;
    };
}
