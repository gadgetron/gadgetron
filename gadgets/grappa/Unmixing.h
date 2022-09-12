#pragma once

#include <map>
#include <memory>

#include "parallel/Merge.h"
#include "Channel.h"

namespace Gadgetron::Grappa {

    class Image {
    public:

        struct {
            uint32_t time_stamp;
            uint16_t slice;
            std::array<float, 3> position, read_dir, phase_dir, slice_dir, table_pos;
        } meta;

        hoNDArray<std::complex<float>> data;
    };

    class Weights {
    public:
        struct {
            uint16_t slice, n_combined_channels, n_uncombined_channels;
        } meta;

        hoNDArray<std::complex<float>> data;
    };

    class Unmixing : public Core::Parallel::Merge {
    public:
        Unmixing(const Core::Context &context, const std::unordered_map<std::string, std::string> &props);

        NODE_PROPERTY(image_series, uint16_t, "Image series number for output images", 0);
        NODE_PROPERTY(unmixing_scale, float, "", 1.0);

        void process(
                std::map<std::string, Core::GenericInputChannel> input,
                Core::OutputChannel output
        ) override;

    private:
        hoNDArray<std::complex<float>> unmix(const Image &image, const Weights &weights);

        std::vector<size_t> create_unmixed_image_dimensions(const Weights &weights);

        static std::vector<size_t> create_output_image_dimensions(const Core::Context &context);
        static std::vector<float> create_output_image_fov(const Core::Context &context);
        ISMRMRD::ImageHeader create_image_header(const Image &image, const Weights &weights);

        const Core::Context context;
        const std::vector<size_t> image_dimensions;
        const std::vector<float> image_fov;

        uint16_t image_index_counter = 0;
    };
}
