#include "ImageArraySplitGadget.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;

namespace {

void doProcessing(AnyImage image, Core::OutputChannel& out) {
    out.push(image);
}

void doProcessing(IsmrmrdImageArray imagearr, Core::OutputChannel& out) {
    
    // 7D, fixed order [X, Y, Z, CHA, N, S, LOC]
    uint16_t X = imagearr.data_.get_size(0);
    uint16_t Y = imagearr.data_.get_size(1);
    uint16_t Z = imagearr.data_.get_size(2);
    uint16_t CHA = imagearr.data_.get_size(3);
    uint16_t N = imagearr.data_.get_size(4);
    uint16_t S = imagearr.data_.get_size(5);
    uint16_t LOC = imagearr.data_.get_size(6);

    // Each image will be [X,Y,Z,CHA] big
    std::vector<size_t> img_dims(4);
    img_dims[0] = X;
    img_dims[1] = Y;
    img_dims[2] = Z;
    img_dims[3] = CHA;

    // Loop over N, S and LOC
    for (uint16_t loc = 0; loc < LOC; loc++) {
        for (uint16_t s = 0; s < S; s++) {
            for (uint16_t n = 0; n < N; n++) {
                // Create a new image header and copy the header for this n, s and loc
                auto imageHeader = ISMRMRD::ImageHeader();
                memcpy(&imageHeader, &imagearr.headers_(n, s, loc), sizeof(ISMRMRD::ImageHeader));

                // Create a new image image
                //  and the 4D data block [X,Y,Z,CHA] for this n, s and loc
                auto imageData = hoNDArray<std::complex<float>>(img_dims);
                memcpy(imageData.get_data_ptr(), &imagearr.data_(0, 0, 0, 0, n, s, loc),
                       X * Y * Z * CHA * sizeof(std::complex<float>));

                // Create a new meta container if needed and copy
                auto imageMetaContainer = std::optional<ISMRMRD::MetaContainer>(); // TODO: Should this be empty? Seems like it should be populated somehow
                if (imagearr.meta_.size() > 0) {
                    size_t mindex = loc * N * S + s * N + n;
                    imageMetaContainer = imagearr.meta_[mindex];
                }

                // Pass the image down the chain
                out.push(Core::Image<std::complex<float>>(std::move(imageHeader), std::move(imageData), std::move(imageMetaContainer)));
            }
        }
    }
}

} // namespace

namespace Gadgetron {

void ImageArraySplitGadget::process(Core::InputChannel<ImageOrImageArray>& in, Core::OutputChannel& out) {

    // Lambda, adds image and correct index to vector of ImageEntries 
    auto lambda = [&](auto message){ 
        doProcessing(message, out);
    };

    for (auto msg : in) {
        visit(lambda, msg);
    }
}

GADGETRON_GADGET_EXPORT(ImageArraySplitGadget);
} // namespace Gadgetron