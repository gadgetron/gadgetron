#pragma once

#include "vector_td.h"
#include <vector>
#include <complex>
#include "hoNDArray.h"
#include <armadillo>
#include <boost/optional.hpp>
#include "mri_core_export.h"

namespace Gadgetron {
    namespace GIRF {


        struct EXPORTMRICORE GIRFdata {

            float sampling_time_us;
            hoNDArray<std::complex<float>> girf_kernel;

        };

        EXPORTMRICORE hoNDArray<std::complex<float>> load_girf_kernel(const std::string& girf_string );


        EXPORTMRICORE hoNDArray<floatd3>
        girf_correct(const hoNDArray<floatd3> &gradients, const hoNDArray <std::complex<float>> &girf_kernel,
                             const arma::fmat33 &rotation_matrix, float gradient_sampling_time,
                             float girf_sampling_time, float TE);

        EXPORTMRICORE hoNDArray<floatd2>
        girf_correct(const hoNDArray<floatd2> &gradients, const hoNDArray <std::complex<float>> &girf_kernel,
                             const arma::fmat33 &rotation_matrix, float gradient_sampling_time,
                             float girf_sampling_time, float TE);
    }
}
