#include "cpu/hoNDArray_fileio.h"
#include "Unmixing.h"

#include "parallel/Merge.h"

#include "hoNDArray.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"


namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Grappa;

    class WeightsProvider {
    public:
        WeightsProvider(
                const Core::Context &context, Core::InputChannel<Weights> &source
        ) : weights(number_of_slices(context), std::nullopt), source(source) {}


        const Weights &operator[](size_t slice) {

            consume_all_pending_weights();
            consume_weights_until_present(slice);

            return weights[slice].value();
        }

    private:

        void consume_all_pending_weights() {
            while (auto w = source.try_pop()) add(w);
        }

        void consume_weights_until_present(size_t slice) {
            while(!weights[slice]) add(source.pop());
        }

        void add(std::optional<Weights> w) {
            weights[w->meta.slice] = std::move(w);
        }

        static size_t number_of_slices(const Core::Context &context) {
            auto e_limits = context.header.encoding[0].encoding_limits;
            return e_limits.slice ? e_limits.slice->maximum + 1u : 1u;
        }

        Core::InputChannel<Weights> &source;
        std::vector<std::optional<Weights>> weights;
    };
}

namespace Gadgetron::Grappa {
    GADGETRON_MERGE_EXPORT(Unmixing);

    Unmixing::Unmixing(
            const Core::Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : Merge(props), context(context),
        image_dimensions(create_output_image_dimensions(context)),
        image_fov(create_output_image_fov(context)) {}

    void Unmixing::process(
            std::map<std::string, Core::GenericInputChannel> input,
            Core::OutputChannel output
    ) {
        Core::InputChannel<Image> images(input.at("images"), output);
        Core::InputChannel<Weights> weights(input.at("weights"), output);

        WeightsProvider weights_provider(context, weights);

        for (auto image : images) {
            auto current_weights = weights_provider[image.meta.slice];
            mrd::Image<std::complex<float>> out;
            out.head = create_image_header(image, current_weights);
            out.data = unmix(image, current_weights);
            output.push(std::move(out));
        }
    }

    hoNDArray<std::complex<float>> Unmixing::unmix(
            const Image &image,
            const Weights &weights
    ) {
        hoNDArray<std::complex<float>> unmixed_image(create_unmixed_image_dimensions(weights));
        clear(unmixed_image);

        auto sets = weights.data.get_number_of_elements() / image.data.get_number_of_elements();
        auto image_elements = unmixed_image.get_number_of_elements() / sets;
        auto coils = weights.data.get_number_of_elements() / (sets * image_elements);

        for (unsigned int s = 0; s < sets; s++) {
            for (unsigned int p = 0; p < image_elements; p++) {
                for (unsigned int c = 0; c < coils; c++) {
                    unmixed_image[s * image_elements + p] +=
                            weights.data[s * image_elements * coils + c * image_elements + p] *
                            image.data[c * image_elements + p] *
                            unmixing_scale;
                }
            }
        }

        return std::move(unmixed_image);
    }

    mrd::ImageHeader Unmixing::create_image_header(const Image &image, const Weights &weights) {

        mrd::ImageHeader header;

        header.slice = image.meta.slice;
        header.acquisition_time_stamp = image.meta.time_stamp;

        std::copy(image_fov.begin(), image_fov.end(), std::begin(header.field_of_view));

        std::copy(image.meta.position.begin(), image.meta.position.end(), std::begin(header.position));
        std::copy(image.meta.read_dir.begin(), image.meta.read_dir.end(), std::begin(header.col_dir));
        std::copy(image.meta.phase_dir.begin(), image.meta.phase_dir.end(), std::begin(header.line_dir));
        std::copy(image.meta.slice_dir.begin(), image.meta.slice_dir.end(), std::begin(header.slice_dir));
        std::copy(image.meta.table_pos.begin(), image.meta.table_pos.end(), std::begin(header.patient_table_position));

        header.image_index = ++image_index_counter;
        header.image_series_index = image_series;

        return header;
    }

    std::vector<size_t> Unmixing::create_output_image_dimensions(const Core::Context &context) {
        auto r_space  = context.header.encoding[0].recon_space;
        return {
                r_space.matrix_size.x,
                r_space.matrix_size.y,
                r_space.matrix_size.z
        };
    }

    std::vector<float> Unmixing::create_output_image_fov(const Core::Context &context) {
        auto r_space  = context.header.encoding[0].recon_space;
        return {
                r_space.field_of_view_mm.x,
                r_space.field_of_view_mm.y,
                r_space.field_of_view_mm.z
        };
    }

    std::vector<size_t> Unmixing::create_unmixed_image_dimensions(const Weights &weights) {
        std::vector<size_t> dimensions = image_dimensions;
        dimensions.push_back(1u + weights.meta.n_uncombined_channels);
        return dimensions;
    }
}