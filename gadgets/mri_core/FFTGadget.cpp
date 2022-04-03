#include "FFTGadget.h"
#include "hoNDArray_math.h"
#include "hoNDFFT.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;

namespace {
    template <class T> Image<T> fft(Image<T>& image) {
        GDEBUG("FFTGadget is not well defined for real-valued images. Doing nothing.");
        return image;
    }

    template <class T> Image<std::complex<T>> fft(Image<std::complex<T>>& image) {
        auto& header = std::get<ISMRMRD::ImageHeader>(image);
        auto& data = std::get<hoNDArray<std::complex<T>>>(image);
		auto &meta = std::get<optional<ISMRMRD::MetaContainer>>(image);
        // Do the FFTs in place
        hoNDFFT<T>::instance()->ifft3c(data);
        return Image<std::complex<T>>(header, data, meta);
    }
} // namespace

namespace Gadgetron {
    AnyImage FFTGadget::process_function(AnyImage image) const {
        return visit([&](auto& image) -> AnyImage { return fft(image); }, image);
    }
GADGETRON_GADGET_EXPORT(FFTGadget);
} // namespace Gadgetron
