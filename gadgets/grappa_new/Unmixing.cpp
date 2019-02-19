#include "Unmixing.h"

#include "parallel/Merge.h"

#include "hoNDArray.h"
#include "hoNDArray_elemwise.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    using namespace Gadgetron::Grappa;

    class WeightsProvider {
    public:
        WeightsProvider(
                const Context &context,
                TypedInputChannel<Weights> &source
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

        TypedInputChannel<Weights> &source;
        std::vector<optional<Weights>> weights;
    };

    void unmix(
            const hoNDArray<std::complex<float>> &weights,
            const hoNDArray<std::complex<float>> &data_in,
            hoNDArray<std::complex<float>> data_out,
            std::complex<float> scale = std::complex<float>(1.0, 0.0)
    ) {
        auto sets = weights.get_number_of_elements() / data_in.get_number_of_elements();
        auto image_elements = data_out.get_number_of_elements() / sets;
        auto coils = weights.get_number_of_elements() / (sets * image_elements);

        clear(data_out);

        for (unsigned int s = 0; s < sets; s++) {
            for (unsigned int p = 0; p < image_elements; p++) {
                for (unsigned int c = 0; c < coils; c++) {
                    data_out[s * image_elements + p] +=
                            weights[s * image_elements * coils + c * image_elements + p] *
                            data_in[c * image_elements + p] * scale;
                }
            }
        }
    }
}

namespace Gadgetron::Grappa {
    GADGETRON_MERGE_EXPORT(Unmixing);

    Unmixing::Unmixing(
            const Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : Merge(props), context(context), image_dimensions(create_output_image_dimensions(context)) {}

    void Unmixing::process(
            std::map<std::string, InputChannel> input,
            OutputChannel output
    ) {
        TypedInputChannel<Image> images(input.at("images"), output);
        TypedInputChannel<Weights> weights(input.at("weights"), output);

        WeightsProvider weights_provider(context, weights);

        for (auto image : images) {
            GINFO_STREAM("Received unmixing job. Excellent.");

            auto current_weights = weights_provider[image.meta.slice];

            hoNDArray<std::complex<float>> unmixed_image(image_dimensions);
            unmix(current_weights.data, image.data, unmixed_image, unmixing_scale);

            GINFO_STREAM("All right! Unmixed image!")
        }
    }

    std::vector<size_t> Unmixing::create_output_image_dimensions(const Core::Context &context) {
        auto r_space  = context.header.encoding[0].reconSpace;
        return {
                r_space.matrixSize.x,
                r_space.matrixSize.y,
                r_space.matrixSize.z
        };
    }
}