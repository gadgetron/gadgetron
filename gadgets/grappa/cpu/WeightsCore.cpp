#include "WeightsCore.h"

namespace Gadgetron::Grappa::CPU {

    const hoNDArray<std::complex<float>> &WeightsCore::estimate_coil_map(const hoNDArray<std::complex<float>> &data) {

        hoNDFFT<float>::instance()->ifft2c(data, buffers.image);
        Gadgetron::coil_map_2d_Inati(buffers.image, buffers.coil_map, coil_map_params.ks, coil_map_params.power);

        return buffers.coil_map;
    }

    hoNDArray<std::complex<float>> WeightsCore::fill_in_uncombined_weights(
            hoNDArray<std::complex<float>> &unmixing_coefficients,
            size_t n_combined_channels
    ) {
        std::vector<hoNDArray<std::complex<float>>> weights = { unmixing_coefficients };

        auto uncombined_coefficients_range = spans(buffers.image_domain_kernel, 3);
        std::for_each(
                uncombined_coefficients_range.begin() + n_combined_channels,
                uncombined_coefficients_range.end(),
                [&](auto w) { weights.push_back(w); }
        );

        return concat(weights);
    }

    hoNDArray<std::complex<float>> WeightsCore::calculate_weights(
            const hoNDArray<std::complex<float>> &data,
            std::array<uint16_t, 4> region_of_support,
            uint16_t acceleration_factor,
            uint16_t n_combined_channels,
            uint16_t n_uncombined_channels
    ) {
        // TODO: Optimize accel_factor == 1;


        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t CHA = data.get_size(2);



        auto coil_map = estimate_coil_map(data);


        Gadgetron::grappa2d_calib_convolution_kernel(
                data,
                data,
                acceleration_factor,
                kernel_params.threshold,
                kernel_params.width,
                kernel_params.height,
                region_of_support[0],
                region_of_support[1],
                region_of_support[2],
                region_of_support[3],
                buffers.convolution_kernel
        );

        Gadgetron::grappa2d_image_domain_kernel(
                buffers.convolution_kernel,
                RO,
                E1,
                buffers.image_domain_kernel
        );

        hoNDArray<std::complex<float>> unmixing_coefficients(RO, E1, CHA);

        Gadgetron::grappa2d_unmixing_coeff(
                buffers.image_domain_kernel,
                coil_map,
                acceleration_factor,
                unmixing_coefficients,
                buffers.g_factor
        );

        return fill_in_uncombined_weights(
                unmixing_coefficients,
                n_combined_channels
        );
    }
}