#pragma once

#include "hoNDArray.h"
#include "denoise_export.h"


namespace Gadgetron {
    namespace Denoise {
        EXPORTDENOISE hoNDArray<float> non_local_means(const hoNDArray<float>& image, float noise_std, unsigned int search_radius);
    }
}
