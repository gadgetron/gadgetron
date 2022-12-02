#pragma once

#include "hoNDArray.h"
#include "cuNDArray.h"
#include "cuNDFFT.h"

namespace Gadgetron::Grappa::GPU {

    class WeightsCore {
    public:
        hoNDArray<std::complex<float>> calculate_weights(
                const hoNDArray<std::complex<float>> &data,
                std::array<uint16_t, 4> region_of_support,
                uint16_t acceleration_factor,
                uint16_t n_combined_channels,
                uint16_t n_uncombined_channels
        );

        cuNDArray<complext<float>> estimate_coil_map(const cuNDArray<complext<float>> &);

        struct {
            uint16_t ks, power;
        } coil_map_params;

        struct {
            uint16_t width, height;
            float threshold;
        } kernel_params;


    };
}