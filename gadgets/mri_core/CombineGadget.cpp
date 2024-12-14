#include "CombineGadget.h"

using namespace Gadgetron;

namespace {
    template <class T> mrd::Image<T> combine(mrd::Image<T>& image) {
        GDEBUG("CombineGadget is not well defined for real-valued images. Doing nothing.");
        return image;
    }

    template <class T>
    mrd::Image<std::complex<T>> combine(mrd::Image<std::complex<T>>& image) {
        // Get the dimensions
        size_t nx = image.data.get_size(0);
        size_t ny = image.data.get_size(1);
        size_t nz = image.data.get_size(2);
        size_t nc = image.data.get_size(3);

        std::vector<size_t> dimensions(4);
        dimensions[0] = nx;
        dimensions[1] = ny;
        dimensions[2] = nz;
        dimensions[3] = 1;

        hoNDArray<std::complex<T>> combinedData(dimensions);

        std::complex<T>* source = image.data.get_data_ptr();
        std::complex<T>* combined = combinedData.get_data_ptr();

        size_t img_block = nx * ny * nz;

        for (size_t z = 0; z < nz; z++) {
            for (size_t y = 0; y < ny; y++) {
                for (size_t x = 0; x < nx; x++) {
                    float mag = 0;
                    float phase = 0;
                    size_t offset = z * ny * nx + y * nx + x;
                    for (size_t c = 0; c < nc; c++) {
                        float mag_tmp = norm(source[offset + c * img_block]);
                        phase += mag_tmp * arg(source[offset + c * img_block]);
                        mag += mag_tmp;
                    }
                    combined[offset] = std::polar(std::sqrt(mag), phase / mag);
                }
            }
        }

        image.data = combinedData;
        return image;
    }
} // namespace

namespace Gadgetron {
    mrd::AnyImage CombineGadget::process_function(mrd::AnyImage image) const {
        return visit([&](auto& image) -> mrd::AnyImage { return combine(image); }, image);
    }
    GADGETRON_GADGET_EXPORT(CombineGadget);
} // namespace Gadgetron
