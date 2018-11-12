#pragma once

#include "vector_td.h"
#include <vector>
#include <complex>
#include "hoNDArray.h"
#include <armadillo>
#include <boost/optional.hpp>

namespace Gadgetron {
    namespace GIRF {


        struct GIRFdata {

            float sampling_time_us;
            hoNDArray<std::complex<float>> girf_kernel;

        };

        hoNDArray<std::complex<float>> load_girf_kernel(const std::string& girf_string );


        hoNDArray<floatd3>
        girf_correct(const hoNDArray<floatd3> &gradients, const hoNDArray <std::complex<float>> &girf_kernel,
                             const arma::fmat33 &rotation_matrix, float gradient_sampling_time,
                             float girf_sampling_time, float TE);

        hoNDArray<floatd2>
        girf_correct(const hoNDArray<floatd2> &gradients, const hoNDArray <std::complex<float>> &girf_kernel,
                             const arma::fmat33 &rotation_matrix, float gradient_sampling_time,
                             float girf_sampling_time, float TE);
    }
}
