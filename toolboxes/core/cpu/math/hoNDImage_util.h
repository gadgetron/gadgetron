/** \file hoNDImage_util.h
\brief math operations on the hoNDImage class.
*/

#pragma once

#include "ho2DArray.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include "ho5DArray.h"
#include "ho6DArray.h"
#include "ho7DArray.h"
#include "hoNDImage.h"
#include "cpucore_math_export.h"

#include "GadgetronCommon.h"
#include <complex>

#include "hoNDArray_reductions.h"
#include "hoNDArray_math_util.h"
#include "hoNDInterpolator.h"

namespace Gadgetron
{
    /// compute the gradient for an ND image
    /// the central difference is computed, the border-value boundary condition is used
    template<class T, unsigned int D> EXPORTCPUCOREMATH bool gradient(const hoNDImage<T, D>& x, hoNDImage<T, D> gx[]);

    /// compute a gaussian kernel
    template<class T> EXPORTCPUCOREMATH bool gaussianKernel(T sigma, double kerWidthInUnitOfSigma, double deltaKer, hoNDArray<T>& ker);

    /// perform the gaussian filter for every dimension
    /// sigma is in the unit of pixel
    template<class ArrayType, class T2> EXPORTCPUCOREMATH bool filterGaussian(ArrayType& x, T2 sigma[], typename ArrayType::value_type* mem=NULL);

    /// perform midian filter
    /// w is the window size
    template<class ArrayType> bool filterMedian(const ArrayType& img, size_t w[], ArrayType& img_out);

    /// downsample the image by a ratio
    /// new image size = image size / ratio
    /// e.g., if ratio = 2, downsample by 2
    template<typename T, typename InterpolatorType, unsigned int D> 
    bool downsampleImage(const hoNDImage<T, D>& in, InterpolatorType& interp, hoNDImage<T, D>& out, float ratio[]);

    /// upsample the image by a ratio
    /// new image size = image size * ratio
    /// e.g., if ratio = 2, upsample by 2
    template<typename T, typename InterpolatorType, unsigned int D> 
    bool upsampleImage(const hoNDImage<T, D>& in, InterpolatorType& interp, hoNDImage<T, D>& out, float ratio[]);

    /// resample the image to specific image size
    /// input and output images occupy the same space region
    /// the pixel size of output images are adjusted accordingly
    template<typename T, typename InterpolatorType, unsigned int D> 
    bool resampleImage(const hoNDImage<T, D>& in, InterpolatorType& interp, const std::vector<size_t>& dim_out, hoNDImage<T, D>& out);

    /// reduce image size by 2 with averaging across two neighbors
    template<typename T, typename BoundaryHandlerType, unsigned int D> 
    bool downsampleImageBy2WithAveraging(const hoNDImage<T, D>& in, BoundaryHandlerType& bh, hoNDImage<T, D>& out);

    /// expand image size by 2 with linear interpolation
    template<typename T, typename BoundaryHandlerType, unsigned int D> 
    bool expandImageBy2(const hoNDImage<T, D>& in, BoundaryHandlerType& bh, hoNDImage<T, D>& out);

    /// filter the image along the first dimension using a 1D kernel
    template<class ArrayType> bool filter1D(const ArrayType& img, const hoNDArray<typename realType<typename ArrayType::value_type>::Type>& ker, GT_BOUNDARY_CONDITION bh, ArrayType& img_out);

    /**
    * @brief r = correlation_coefficient(a, b)
    */
    template <typename T, unsigned int D> 
    bool corrCoef(const hoNDImage<T, D>& a, const hoNDImage<T, D>& b, T& r);
}

#include "hoNDImage_util.hxx"
