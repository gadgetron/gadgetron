#pragma once

#include <map>
#include <memory>

#include "parallel/Merge.h"
#include "Channel.h"

namespace Gadgetron::Grappa {

    class Image {
    public:

        struct {
            uint64_t slice;
            std::array<float, 3> position, read_dir, phase_dir, slice_dir, table_pos;
        } meta;

        hoNDArray<std::complex<float>> data;
    };

    class Weights {
    public:
        struct {
            uint64_t slice, n_combined_channels;
        } meta;

        hoNDArray<std::complex<float>> data;
    };

    class Unmixing : public Core::Parallel::Merge {
    public:
        Unmixing(const Core::Context &context, const std::unordered_map<std::string, std::string> &props);


        // TODO: Get rid of these pointless things.
        NODE_PROPERTY(device_channels, int, "Number of device channels", 0);

        NODE_PROPERTY(image_series, int, "Image series number for output images", 0);

        NODE_PROPERTY(target_channels, int, "Number of target channels for GRAPPA recon", 0);


        void process(
                std::map<std::string, Core::InputChannel> input,
                Core::OutputChannel output
        ) override;
    };
}
