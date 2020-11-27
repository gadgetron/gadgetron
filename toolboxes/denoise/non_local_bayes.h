//
// Created by dchansen on 6/20/18.
//
#pragma once

#include "hoNDArray.h"

namespace Gadgetron {
    namespace Denoise {
        hoNDArray<float> non_local_bayes(const hoNDArray<float>& image, float noise_std=1.0f, unsigned int search_radius=25);
        hoNDArray<std::complex<float>> non_local_bayes(const hoNDArray<std::complex<float>>& image, float noise_std=1.0f, unsigned int search_radius=25);
    }
}
