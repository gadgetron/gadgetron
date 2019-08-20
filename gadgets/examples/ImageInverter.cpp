#include "ImageInverter.h"

#include "hoNDArray_math.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;

namespace {
    template<class T>
    Image<T> invert_image(const Image<T> &image) {

        auto data = std::get<hoNDArray<T>>(image);

        data -= max(data);
        data *= T(-1);

        return Image<T>(
                std::get<ISMRMRD::ImageHeader>(image),
                data,
                std::get<optional<ISMRMRD::MetaContainer>>(image)
        );
    }

    template<class T>
    Image<std::complex<T>> invert_image(const Image<std::complex<T>> &image) {
        throw std::runtime_error("Invert image is not well defined for complex images.");
    }
}

namespace Gadgetron::Examples {

    ImageInverter::ImageInverter(const Context &, const GadgetProperties &properties) : TypedPureGadget(properties) {}

    AnyImage ImageInverter::process_function(AnyImage image) const {
        GINFO_STREAM("Inverting image.")
        return apply_visitor([](const auto &image) -> AnyImage { return invert_image(image); }, image);
    }

    GADGETRON_GADGET_EXPORT(ImageInverter);
}

