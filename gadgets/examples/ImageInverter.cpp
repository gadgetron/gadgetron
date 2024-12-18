#include "ImageInverter.h"

#include "hoNDArray_math.h"

using namespace Gadgetron;

namespace {
    template<class T>
    mrd::Image<T> invert_image(const mrd::Image<T> &image)
    {
        mrd::Image<T> out;
        out.head = image.head;
        out.meta = image.meta;
        out.data = image.data;

		auto max_value = *std::max_element(image.data.begin(), image.data.end());
		for (auto& d : out.data) d = max_value - d;

        return out;
    }

    template<class T>
    mrd::Image<std::complex<T>> invert_image(const mrd::Image<std::complex<T>> &image) {
        throw std::runtime_error("Invert image is not well defined for complex images.");
    }
}

namespace Gadgetron::Examples {

    mrd::AnyImage ImageInverter::process_function(mrd::AnyImage image) const {
        GINFO_STREAM("Inverting image.")
        return visit([](const auto &image) -> mrd::AnyImage { return invert_image(image); }, image);
    }

    GADGETRON_GADGET_EXPORT(ImageInverter);
}

