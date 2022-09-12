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

#include <complex>

#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNDInterpolator.h"

namespace Gadgetron
{
    /// compute the gradient for an ND image
    /// the central difference is computed, the border-value boundary condition is used
    template <typename ImageType> bool gradient(const ImageType& x, ImageType gx[]);

    /// compute a gaussian kernel
    template<class T> bool gaussianKernel(T sigma, double kerWidthInUnitOfSigma, double deltaKer, hoNDArray<T>& ker);

    /// perform the gaussian filter for every dimension
    /// sigma is in the unit of pixel
    template<class ArrayType, class T2> bool filterGaussian(ArrayType& x, T2 sigma[], typename ArrayType::value_type* mem=NULL);

    /// perform midian filter
    /// w is the window size
    template<class ArrayType> bool filterMedian(const ArrayType& img, size_t w[], ArrayType& img_out);

    /// downsample the image by a ratio
    /// new image size = image size / ratio
    /// e.g., if ratio = 2, downsample by 2
    template<typename ImageType, typename InterpolatorType> 
    bool downsampleImage(const ImageType& in, InterpolatorType& interp, ImageType& out, float ratio[]);

    /// upsample the image by a ratio
    /// new image size = image size * ratio
    /// e.g., if ratio = 2, upsample by 2
    template<typename ImageType, typename InterpolatorType> 
    bool upsampleImage(const ImageType& in, InterpolatorType& interp, ImageType& out, float ratio[]);

    /// resample the image to specific image size
    /// input and output images occupy the same space region
    /// the pixel size of output images are adjusted accordingly
    template<typename ImageType, typename InterpolatorType> 
    bool resampleImage(const ImageType& in, InterpolatorType& interp, const std::vector<size_t>& dim_out, ImageType& out);

    /// reduce image size by 2 with averaging across two neighbors
    template<typename ImageType, typename BoundaryHandlerType> 
    bool downsampleImageBy2WithAveraging(const ImageType& in, BoundaryHandlerType& bh, ImageType& out);

    /// expand image size by 2 with linear interpolation
    template<typename ImageType, typename BoundaryHandlerType> 
    bool expandImageBy2(const ImageType& in, BoundaryHandlerType& bh, ImageType& out);

    /// filter the image along the first dimension using a 1D kernel
    template<class ArrayType> bool filter1D(const ArrayType& img, const hoNDArray<typename realType<typename ArrayType::value_type>::Type>& ker, GT_BOUNDARY_CONDITION bh, ArrayType& img_out);

    /**
    * @brief r = correlation_coefficient(a, b)
    */
    template <typename ImageType> 
    bool corrCoef(const ImageType& a, const ImageType& b, typename ImageType::value_type& r);
}

#include "hoNDImage_util.hxx"
