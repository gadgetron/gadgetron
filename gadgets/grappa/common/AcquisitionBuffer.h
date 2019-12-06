#pragma once

#include <map>
#include <set>
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
        view(size_t index) const;

        void clear(size_t index);

        bool is_fully_sampled(size_t index) const;

        void add_pre_update_callback(std::function<void(const AnnotatedAcquisition &)> fn);
        void add_post_update_callback(std::function<void(const AnnotatedAcquisition &)> fn);

    private:
        const Core::Context context;

        struct {
            std::set<uint32_t> expected_lines;
            size_t number_of_samples, number_of_lines;
            int line_offset = 0;
        } internals;

        struct buffer {
            hoNDArray<std::complex<float>> data;
            std::set<uint32_t> sampled_lines;
        };

        buffer create_buffer(const std::vector<size_t> &dimensions);

        std::map<size_t, buffer> buffers;

        std::vector<std::function<void(const AnnotatedAcquisition &)>> pre_update_callbacks;
        std::vector<std::function<void(const AnnotatedAcquisition &)>> post_update_callbacks;
    };
}