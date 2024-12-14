#include "ImageArraySplitGadget.h"


using namespace Gadgetron;

namespace {

void splitInputData(mrd::AnyImage image, Core::OutputChannel& out) {
    out.push(image);
}

void splitInputData(mrd::ImageArray imagearr, Core::OutputChannel& out) {
        // 7D, fixed order [X, Y, Z, CHA, N, S, LOC]
    uint16_t X = imagearr.data.get_size(0);
    uint16_t Y = imagearr.data.get_size(1);
    uint16_t Z = imagearr.data.get_size(2);
    uint16_t CHA = imagearr.data.get_size(3);
    uint16_t N = imagearr.data.get_size(4);
    uint16_t S = imagearr.data.get_size(5);
    uint16_t LOC = imagearr.data.get_size(6);

    // Each image will be [X,Y,Z,CHA] big
    std::vector<size_t> img_dims(4);
    img_dims[0] = X;
    img_dims[1] = Y;
    img_dims[2] = Z;
    img_dims[3] = CHA;

    // Loop over LOC, S, and N
    for (auto loc = 0; loc < LOC; loc++) {
        for (auto s = 0; s < S; s++) {
            for (auto n = 0; n < N; n++) {
                mrd::Image<std::complex<float>> img;
                img.head = imagearr.headers(n, s, loc);
                if (imagearr.meta.size() >= LOC * S * N) {
                    img.meta = imagearr.meta(n, s, loc);
                }

                img.data.create(img_dims);

                memcpy(img.data.data(), &imagearr.data(0, 0, 0, 0, n, s, loc), X * Y * Z * CHA * sizeof(std::complex<float>));

                // Pass the image down the chain
                out.push(std::move(img));
            }
        }
    }
}

} // namespace

namespace Gadgetron {

void ImageArraySplitGadget::process(Core::InputChannel<ImageOrImageArray>& in, Core::OutputChannel& out) {
    for (auto msg : in) {
        visit([&](auto message){splitInputData(message, out);}, msg);
    }
}

GADGETRON_GADGET_EXPORT(ImageArraySplitGadget);
} // namespace Gadgetron
