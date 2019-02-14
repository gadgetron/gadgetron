#pragma once

#include <functional>

#include "AnnotatedAcquisition.h"

#include "Context.h"
#include "Channel.h"
#include "Types.h"

#include "hoNDArray.h"

namespace Gadgetron::Grappa {

    class AcquisitionBuffer {
    public:
        explicit AcquisitionBuffer(Core::Context);

        void add(const AnnotatedAcquisition &acquisition);

        template<class T>
        void add(const T &acquisitions) {
            for (const auto &acquisition : acquisitions) {
                add(acquisition);
            }
        }

        hoNDArray<std::complex<float>>
        take(size_t index);

        const hoNDArray<std::complex<float>> &
        view(size_t index);

        void clear(size_t index);

        void add_pre_update_callback(std::function<void(const AnnotatedAcquisition &)> fn);
        void add_post_update_callback(std::function<void(const AnnotatedAcquisition &)> fn);

    private:
        const Core::Context context;

        struct {
            std::vector<unsigned long> buffer_dimensions;
            int line_offset = 0;
        } internals;

        std::vector<hoNDArray<std::complex<float>>> buffers;

        std::vector<std::function<void(const AnnotatedAcquisition &)>> pre_update_callbacks;
        std::vector<std::function<void(const AnnotatedAcquisition &)>> post_update_callbacks;
    };
}