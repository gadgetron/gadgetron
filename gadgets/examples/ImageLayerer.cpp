#include "ImageLayerer.h"

#include <vector>

#include "hoNDArray_math.h"
#include "hoNDArray_utils.h"

#include "mri_core_utility.h"

using namespace Gadgetron;
namespace {

    template<class T>
    mrd::Image<T> merge(const mrd::Image<T> &a, const mrd::Image<T> &b)
    {
        const auto& data_a = a.data;
        const auto& data_b = b.data;

        if (data_a.dimensions() != data_b.dimensions()) {
            throw std::runtime_error("Images missized. Can't merge mismatched images.");
        }

        std::vector<size_t> size = {
                data_a.get_size(0),
                data_a.get_size(1),
                data_a.get_size(2),
                data_a.get_size(4) * 2
        };

        auto data = concat(std::vector<hoNDArray<T>>{data_a, data_b});
        data.reshape(size);

        mrd::Image<T> out{.head=a.head, .data=data, .meta=a.meta};

        return out;
    }

    template<class A, class B>
    mrd::Image<A> merge(const mrd::Image<A> &, const mrd::Image<B> &) {
        throw std::runtime_error("Images have different types; merging will be a hassle.");
    }
}

namespace Gadgetron::Examples {

    ImageLayerer::ImageLayerer(const Core::Context &, const Core::GadgetProperties &properties) : Merge(properties) {}

    void ImageLayerer::process(std::map<std::string, Core::GenericInputChannel> input, Core::OutputChannel output) {

        auto unchanged = Core::InputChannel<mrd::AnyImage>(input.at("unchanged"), output);
        auto inverted = Core::InputChannel<mrd::AnyImage>(input.at("inverted"), output);

        for (auto image : unchanged) {
            auto merged = std::visit(
                    [](const auto &a, const auto &b) -> mrd::AnyImage { return merge(a, b); },
                    image,
                    inverted.pop()
            );



            GINFO_STREAM("Images combined; pushing out result.");

            output.push(std::move(merged));
        }
    }

    GADGETRON_MERGE_EXPORT(ImageLayerer)
}
