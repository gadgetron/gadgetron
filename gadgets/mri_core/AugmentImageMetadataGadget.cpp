#include "AugmentImageMetadataGadget.h"

#include "log.h"

namespace Gadgetron{

mrd::Image<std::complex<float>> AugmentImageMetadataGadget::process_function(
    mrd::Image<std::complex<float>> input_image) const
{
    if (input_image.meta["ImageRowDir"].size() != 3) {
        input_image.meta["ImageRowDir"].push_back(input_image.head.col_dir[0]);
        input_image.meta["ImageRowDir"].push_back(input_image.head.col_dir[1]);
        input_image.meta["ImageRowDir"].push_back(input_image.head.col_dir[2]);
    }


    if (input_image.meta["ImageColumnDir"].size() != 3) {
        input_image.meta["ImageColumnDir"].push_back(input_image.head.line_dir[0]);
        input_image.meta["ImageColumnDir"].push_back(input_image.head.line_dir[1]);
        input_image.meta["ImageColumnDir"].push_back(input_image.head.line_dir[2]);
    }

    return input_image;
}

GADGETRON_GADGET_EXPORT(AugmentImageMetadataGadget)
}