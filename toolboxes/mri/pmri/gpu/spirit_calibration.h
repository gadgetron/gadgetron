/** \file spirit_calibration.h
    \brief Utility to calibrate spirit convolution kernels, GPU-based.
*/

#pragma once

#include "cuNDArray.h"

namespace Gadgetron
{

  /**
     @brief Utility to estimate spirit convolution kernels, GPU-based.
     @param[in] cartesian_kspace_data Array with fully sampled kspace data (Cartesian). E.g. as a result of accumulation of multiple frames.
     @param[in] kernel_size Size of the convolution kernel to use for k-space calibration. Must be an odd number.
     @return A set convolution kernels Fourier transformed into image space. For 'n' coils, n^2 calibration images are estimated, i.e. 'n' kernels for each coil.
     Currently only 2D Spirit is supported in this function (higher-dimensional Spirit is supported in the gt-plus toolbox).
  */
  boost::shared_ptr< cuNDArray<float_complext> >
  estimate_spirit_kernels( cuNDArray<float_complext> *cartesian_kspace_data, unsigned int kernel_size );
}
