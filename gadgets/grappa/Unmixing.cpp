#include "cpu/hoNDArray_fileio.h"
#include "Unmixing.h"

#include "parallel/Merge.h"

#include "hoNDArray.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"


namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Grappa;

    class WeightsProvider {
    public:
        WeightsProvider(
                const Context &context, InputChannel<Weights> &source
        ) : weights(number_of_slices(context), none), source(source) {}


        const Weights &operator[](size_t slice) {

            consume_all_pending_weights();
            consume_weights_until_present(slice);

            return weights[slice].get();
        }

    private:

        void consume_all_pending_weights() {
            while (auto w = source.try_pop()) add(w);
        }

        void consume_weights_until_present(size_t slice) {
            while(!weights[slice]) add(source.pop());
        }

        void add(optional<Weights> w) {
            weights[w->meta.slice] = std::move(w);
        }

        static size_t number_of_slices(const Context &context) {
            auto e_limits = context.header.encoding[0].encodingLimits;
            return e_limits.slice ? e_limits.slice->maximum + 1u : 1u;
        }

        InputChannel<Weights> &source;
        std::vector<optional<Weights>> weights;
    };
}

namespace Gadgetron::Grappa {
    GADGETRON_MERGE_EXPORT(Unmixing);

    Unmixing::Unmixing(
            const Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : Merge(props), context(context),
        image_dimensions(create_output_image_dimensions(context)),
        image_fov(create_output_image_fov(context)) {}

    void Unmixing::process(
            std::map<std::string, GenericInputChannel> input,
            OutputChannel output
    ) {
        InputChannel<Image> images(input.at("images"), output);
        InputChannel<Weights> weights(input.at("weights"), output);

        WeightsProvider weights_provider(context, weights);

        for (auto image : images) {
            auto current_weights = weights_provider[image.meta.slice];
            output.push(
                    create_image_header(image, current_weights),
                    unmix(image, current_weights)
            );
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

    ISMRMRD::ImageHeader Unmixing::create_image_header(const Image &image, const Weights &weights) {

        ISMRMRD::ImageHeader header;

        header.slice = image.meta.slice;
        header.acquisition_time_stamp = image.meta.time_stamp;
        header.channels = weights.meta.n_uncombined_channels + uint16_t(1);

        std::copy(image_dimensions.begin(), image_dimensions.end(), std::begin(header.matrix_size));
        std::copy(image_fov.begin(), image_fov.end(), std::begin(header.field_of_view));

        std::copy(image.meta.position.begin(), image.meta.position.end(), std::begin(header.position));
        std::copy(image.meta.read_dir.begin(), image.meta.read_dir.end(), std::begin(header.read_dir));
        std::copy(image.meta.phase_dir.begin(), image.meta.phase_dir.end(), std::begin(header.phase_dir));
        std::copy(image.meta.slice_dir.begin(), image.meta.slice_dir.end(), std::begin(header.slice_dir));
        std::copy(image.meta.table_pos.begin(), image.meta.table_pos.end(), std::begin(header.patient_table_position));

        header.image_index = ++image_index_counter;
        header.image_series_index = image_series;

        return header;
    }

    std::vector<size_t> Unmixing::create_output_image_dimensions(const Core::Context &context) {
        auto r_space  = context.header.encoding[0].reconSpace;
        return {
                r_space.matrixSize.x,
                r_space.matrixSize.y,
                r_space.matrixSize.z
        };
    }

    std::vector<float> Unmixing::create_output_image_fov(const Core::Context &context) {
        auto r_space  = context.header.encoding[0].reconSpace;
        return {
                r_space.fieldOfView_mm.x,
                r_space.fieldOfView_mm.y,
                r_space.fieldOfView_mm.z
        };
    }

    std::vector<size_t> Unmixing::create_unmixed_image_dimensions(const Weights &weights) {
        std::vector<size_t> dimensions = image_dimensions;
        dimensions.push_back(1u + weights.meta.n_uncombined_channels);
        return dimensions;
    }
}