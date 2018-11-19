#pragma once

#include "mri_core_data.h"
namespace Gadgetron { namespace NonCartesian {

        void append_image_header(IsmrmrdImageArray &res, const IsmrmrdReconBit &recon_bit,
                                                          size_t encoding);
    }
}