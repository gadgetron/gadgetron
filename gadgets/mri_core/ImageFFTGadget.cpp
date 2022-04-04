#include "ImageFFTGadget.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;

namespace {
    template<class T>
    Image<std::complex<T>> imageFFT(Image<std::complex<T>> &image) {
        ISMRMRD::ImageHeader &header = std::get<ISMRMRD::ImageHeader>(image);
        hoNDArray<std::complex<T>> &data = std::get<hoNDArray<std::complex<T>>>(image);
		optional<ISMRMRD::MetaContainer> &meta = std::get<optional<ISMRMRD::MetaContainer>>(image);
		
        std::vector<size_t> dims;
        data.get_dimensions(dims);

        size_t RO = dims[0];
        size_t E1 = dims[1];
        size_t E2 = 1;

        if (E2 > 1)
        {
            Gadgetron::hoNDFFT<T>::instance()->fft3c(data);
        }
        else
        {
            Gadgetron::hoNDFFT<T>::instance()->fft2c(data);
        }
        if(meta.has_value()){
            meta.value().append(GADGETRON_IMAGEPROCESSINGHISTORY, "FFT");
        }
        return image;
    }
}

namespace Gadgetron {

    ComplexImage ImageFFTGadget::process_function(ComplexImage image) const {
        return visit([&](auto &image) -> ComplexImage { return imageFFT(image); }, image);
    }

    GADGETRON_GADGET_EXPORT(ImageFFTGadget);
}

