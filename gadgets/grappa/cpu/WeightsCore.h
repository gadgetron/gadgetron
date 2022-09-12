#pragma once

#include "hoNDFFT.h"
#include "hoNDArray.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_iterators.h"
#include "mri_core_grappa.h"
#include "mri_core_coil_map_estimation.h"

namespace Gadgetron::Grappa::CPU {

    class WeightsCore {
    public:
        hoNDArray<std::complex<float>> calculate_weights(
                const hoNDArray<std::complex<float>> &data,
                std::array<uint16_t, 4> region_of_support,
                uint16_t acceleration_factor,
                uint16_t n_combined_channels,
                uint16_t n_uncombined_channels
        );

        const hoNDArray<std::complex<float>> &
        estimate_coil_map(
                const hoNDArray<std::complex<float>> &data
        );

        hoNDArray<std::complex<float>>
        fill_in_uncombined_weights(
                hoNDArray<std::complex<float>> &unmixing_coefficients,
                size_t n_combined_channels
        );

        struct {
            uint16_t ks, power;
        } coil_map_params;

        struct {
            uint16_t width, height;
            float threshold;
        } kernel_params;

        struct {
            // We maintain a few buffers to avoid reallocating them repeatedly.
            hoNDArray<std::complex<float>> image, coil_map, convolution_kernel, image_domain_kernel;
            hoNDArray<float> g_factor;
        } buffers;
    };
}