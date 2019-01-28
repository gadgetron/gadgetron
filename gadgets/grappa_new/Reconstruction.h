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

        NODE_PROPERTY(image_series, int, "Image series number for output images", 0);

        // TODO: This doesn't need to be special? We can handle all of it in 'uncombined_channels'?
        NODE_PROPERTY(uncombined_channels, std::string, "Uncombined channels (as a comma separated list of channel indices", "");
        NODE_PROPERTY(uncombined_channels_by_name, std::string, "Uncombined channels (as a comma separated list of channel names", "");

        void process(
                std::map<std::string, std::shared_ptr<Core::Channel>> input,
                std::shared_ptr<Core::Channel> output
        ) override;
    };
}
