#include "ImageLayerer.h"

#include <vector>

#include "hoNDArray_math.h"
#include "hoNDArray_utils.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;

namespace {

    template<class T>
    Image<T> merge(const Image<T> &a, const Image<T> &b) {

        auto header = std::get<ISMRMRD::ImageHeader>(a);
        const auto &data_a = std::get<hoNDArray<T>>(a);
        const auto &data_b = std::get<hoNDArray<T>>(b);

        header.channels *= 2;

        if (data_a.dimensions() != data_b.dimensions()) {
            throw std::runtime_error("Images missized. Can't merge mismatched images.");
        }

        std::vector<size_t> size = {
                data_a.get_size(0),
                data_a.get_size(1),
                data_a.get_size(2),
                header.channels
        };

        auto data = concat(std::vector<hoNDArray<T>>{data_a, data_b});
        data.reshape(size);

//        auto data = hoNDArray<T>(size);
//        data(slice, slice, 0) = data_a;
//        data(slice, slice, 1) = data_b;

        return Image<T>(header, data, Core::none);
    }

    template<class A, class B>
    Image<A> merge(const Image<A> &, const Image<B> &) {
        throw std::runtime_error("Images have different types; merging will be a hazzle.");
    }
}

namespace Gadgetron::Examples {

    ImageLayerer::ImageLayerer(const Context &, const GadgetProperties &properties) : Merge(properties) {}

    void ImageLayerer::process(std::map<std::string, GenericInputChannel> input, OutputChannel output) {

        auto unchanged = InputChannel<AnyImage>(input.at("unchanged"), output);
        auto inverted = InputChannel<AnyImage>(input.at("inverted"), output);

        for (auto image : unchanged) {
            auto merged = Core::visit(
                    [](const auto &a, const auto &b) -> AnyImage { return merge(a, b); },
                    image,
                    inverted.pop()
            );



            GINFO_STREAM("Images combined; pushing out result.");

            output.push(std::move(merged));
        }
    }

    GADGETRON_MERGE_EXPORT(ImageLayerer)
}
