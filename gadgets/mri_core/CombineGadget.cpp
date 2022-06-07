#include "CombineGadget.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;

namespace {
    template <class T> Image<T> combine(Image<T>& image) {
        GDEBUG("CombineGadget is not well defined for real-valued images. Doing nothing.");
        return image;
    }

    template <class T>
    Image<std::complex<T>> combine(Image<std::complex<T>>& image) {
        auto& header = std::get<ISMRMRD::ImageHeader>(image);
        auto& data = std::get<hoNDArray<std::complex<T>>>(image);

        // Get the dimensions
        size_t nx = data.get_size(0);
        size_t ny = data.get_size(1);
        size_t nz = data.get_size(2);
        size_t nc = data.get_size(3);

        // Create a new message with an hoNDArray for the combined image
        hoNDArray<std::complex<T>> m3 = hoNDArray<std::complex<T>>();

        std::vector<size_t> dimensions(3);
        dimensions[0] = nx;
        dimensions[1] = ny;
        dimensions[2] = nz;

        try {
            m3.create(dimensions);
        } catch (std::runtime_error& err) {
            GEXCEPTION(err, "CombineGadget, failed to allocate new array\n");
        }

        std::complex<T>* d1 = data.get_data_ptr();
        std::complex<T>* d2 = m3.get_data_ptr();

        size_t img_block = nx * ny * nz;

        for (size_t z = 0; z < nz; z++) {
            for (size_t y = 0; y < ny; y++) {
                for (size_t x = 0; x < nx; x++) {
                    float mag = 0;
                    float phase = 0;
                    size_t offset = z * ny * nx + y * nx + x;
                    for (size_t c = 0; c < nc; c++) {
                        float mag_tmp = norm(d1[offset + c * img_block]);
                        phase += mag_tmp * arg(d1[offset + c * img_block]);
                        mag += mag_tmp;
                    }
                    d2[offset] = std::polar(std::sqrt(mag), phase / mag);
                }
            }
        }

        // Modify header to match the size and change the type to real
        header.channels = 1;
        return Image<std::complex<T>>(std::get<ISMRMRD::ImageHeader>(image), std::move(m3), std::get<optional<ISMRMRD::MetaContainer>>(image));
    }
} // namespace

namespace Gadgetron {
    AnyImage CombineGadget::process_function(AnyImage image) const {
        return visit([&](auto& image) -> AnyImage { return combine(image); }, image);
    }
    GADGETRON_GADGET_EXPORT(CombineGadget);
} // namespace Gadgetron
