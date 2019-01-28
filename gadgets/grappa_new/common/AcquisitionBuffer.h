#pragma once

#include <functional>

#include "Context.h"
#include "Channel.h"
#include "Types.h"

#include "hoNDArray.h"

namespace Gadgetron::Grappa {

    class AcquisitionBuffer {
    public:
        explicit AcquisitionBuffer(Core::Context);

        void add(const Core::Acquisition &acquisition);

        template<class T>
        void add(const T &acquisitions) {
            for (const auto &acquisition : acquisitions) {
                add(acquisition);
            }
        }

        void add_acquisition_hook(std::function<void(const Core::Acquisition &)> fn);

        hoNDArray<std::complex<float>>
        take(size_t index);

        const hoNDArray<std::complex<float>> &
        view(size_t index);

        void clear(size_t index);

    private:
        const Core::Context context;

        struct {
            std::vector<unsigned long> buffer_dimensions;
            int line_offset = 0;
        } internals;

        std::vector<hoNDArray<std::complex<float>>> buffers;
        std::vector<std::function<void(const Core::Acquisition &)>> acquisition_hooks;
    };
}