#pragma once

#include "vector_td.h"
#include <vector>
#include <complex>
#include "hoNDArray.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_math.h"
#include <boost/optional.hpp>
#include "hoNDArray_fileio.h"
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include <stdio.h>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <boost/filesystem/fstream.hpp>
#include "armadillo"

using namespace Gadgetron;

namespace Gadgetron
{
    namespace GIRF
    {

        hoNDArray<std::complex<float>> load_girf_kernel(const std::string &girf_string);

        hoNDArray<floatd3>
        girf_correct(const hoNDArray<floatd3> &gradients, const hoNDArray<std::complex<float>> &girf_kernel,
                     const arma::fmat33 &rotation_matrix, float gradient_sampling_time,
                     float girf_sampling_time, float TE);

        hoNDArray<floatd2>
        girf_correct(const hoNDArray<floatd2> &gradients, const hoNDArray<std::complex<float>> &girf_kernel,
                     const arma::fmat33 &rotation_matrix, float gradient_sampling_time,
                     float girf_sampling_time, float TE);

        hoNDArray<std::complex<float>> zeropadding(hoNDArray<std::complex<float>> input, int zpadFactor);
        hoNDArray<std::complex<float>> readGIRFKernel(std::string folder);

    } // namespace GIRF
} // namespace nhlbi_toolbox
