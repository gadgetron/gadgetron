//
// Created by dchansen on 6/20/18.
//
#pragma once

#include "hoNDArray.h"
#include "denoise_export.h"


namespace Gadgetron {
    namespace Denoise {
        EXPORTDENOISE hoNDArray<float> non_local_bayes(const hoNDArray<float>& image, float noise_std, unsigned int search_radius);
        EXPORTDENOISE hoNDArray<std::complex<float>> non_local_bayes(const hoNDArray<std::complex<float>>& image, float noise_std, unsigned int search_radius);
    }
}
