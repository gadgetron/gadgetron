#pragma once

#include <map>
#include <memory>

#include "parallel/Merge.h"
#include "Channel.h"

namespace Gadgetron::Grappa {

    class Image {
    public:

        struct {
            size_t slice;
            std::array<float, 3> position, read_dir, phase_dir, slice_dir, table_pos;
        } meta;

        hoNDArray<std::complex<float>> data;
    };

    class Weights {
    public:

    };

    class Reconstruction : public Core::Parallel::Merge {
    public:
        Reconstruction(const Core::Context &context, const std::unordered_map<std::string, std::string> &props);

        void process(
                std::map<std::string, std::shared_ptr<Core::Channel>> input,
                std::shared_ptr<Core::Channel> output
        ) override;
    };
}
