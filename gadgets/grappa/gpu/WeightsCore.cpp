#include "WeightsCore.h"

#include "../common/grappa_common.h"

#include "cuNDArray.h"
#include "cuNDFFT.h"
#include "htgrappa.h"
#include "b1_map.h"

namespace Gadgetron::Grappa::GPU {

    hoNDArray<std::complex<float>> WeightsCore::calculate_weights(
            const hoNDArray<std::complex<float>> &data,
            std::array<uint16_t, 4> region_of_support,
            uint16_t acceleration_factor,
            uint16_t n_combined_channels,
            uint16_t n_uncombined_channels
    ) {
        if (n_uncombined_channels) {
            throw std::runtime_error("GPU RT Grappa does not currently support uncombined channels.");
        }

        // First line is a crime. Look at this. Look at it! Enjoy the crime! TODO: Fight crime!
        cuNDArray<complext<float>> k_space_data(reinterpret_cast<const hoNDArray<complext<float>> &>(data));
        cuNDArray<complext<float>> coil_map = estimate_coil_map(k_space_data);

        cuNDArray<complext<float>> unmixing_coefficients(data.dimensions());
        std::vector<unsigned int> kernel_size = { kernel_params.width, kernel_params.height };
        std::vector<std::pair<unsigned int, unsigned int>> ros = {
                std::make_pair(region_of_support[0], region_of_support[1]),
                std::make_pair(region_of_support[2], region_of_support[3])
        };

        htgrappa_calculate_grappa_unmixing(
                &k_space_data,
                &coil_map,
                acceleration_factor,
                &kernel_size,
                &unmixing_coefficients,
                &ros
        );

        auto output_ptr = unmixing_coefficients.to_host();
        return std::move(reinterpret_cast<hoNDArray<std::complex<float>>&>(*output_ptr));
    }

    cuNDArray<complext<float>> WeightsCore::estimate_coil_map(const cuNDArray<complext<float>> &k_space_data) {

        cuNDArray<complext<float>> r_space_data(k_space_data);

        std::vector<size_t> fft_dims = {0, 1};
        cuNDFFT<float>::instance()->ifft(&r_space_data, &fft_dims);

        return estimate_b1_map<float, 2>(r_space_data);
    }
}