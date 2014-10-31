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
#include "hoNDMath_util.h"

#include "GadgetronCommon.h"
#include <complex>

#include "hoNDArray_math_util.h"
#include "hoNDInterpolator.h"

namespace Gadgetron
{

    /**
    * @brief Construct a complex array from a real array.
    * @param[in] x Input array.
    * @return A new complex array containing the input array in the real component and zeros in the imaginary component.
    */

    template<class T, unsigned int D> 
    bool real_imag_to_complex(const hoNDImage<typename realType<T>::Type, D>& real, 
                            const hoNDImage<typename realType<T>::Type, D>& imag, 
                            hoNDImage<T, D>& cplx);

    template<class T, unsigned int D> 
    bool complex_to_real_imag(const hoNDImage<T, D>& cplx, 
                            hoNDImage<typename realType<T>::Type, D>& real, 
                            hoNDImage<typename realType<T>::Type, D>& imag);

    template<unsigned int D> 
    bool complex_to_real_imag(const hoNDImage<float, D>& cplx, 
                            hoNDImage<float, D>& real, 
                            hoNDImage<float, D>& imag);

    template<unsigned int D> 
    bool complex_to_real_imag(const hoNDImage<double, D>& cplx, 
                            hoNDImage<double, D>& real, 
                            hoNDImage<double, D>& imag);

    template<class T, unsigned int D> 
    bool complex_to_real(const hoNDImage<T, D>& cplx, hoNDImage<T, D>& real);

    template<class T, unsigned int D> 
    bool complex_to_real(const hoNDImage<T, D>& cplx, hoNDImage<typename realType<T>::Type, D>& real);

    template<class T, unsigned int D> 
    bool complex_to_real(hoNDImage<T, D>& cplx);

    template<class T, unsigned int D> 
    bool complex_to_imag(const hoNDImage<T, D>& cplx, hoNDImage<typename realType<T>::Type, D>& imag);

    template<class T, unsigned int D> 
    bool real_to_complex(const hoNDImage<typename realType<T>::Type, D>& real, hoNDImage<T, D>& cplx);

    template <class T, unsigned int D>
    bool minValue(const hoNDImage<T, D>& im, T& v);

    template <class T, unsigned int D>
    bool maxValue(const hoNDImage<T, D>& im, T& v);

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
    * @brief add two vectors of values, r = x + y
    support in-place computation, e.g. x==r or y==r
    */
    template <typename T, unsigned int D> 
    bool add(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r);

    template <typename T, unsigned int D> 
    bool add(size_t N, const T* x, const T* y, T* r);

    /**
    * @brief subtract two vectors of values, r = x - y
    support in-place computation, e.g. x==r
    */
    template <typename T, unsigned int D> 
    bool subtract(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r);

    template <typename T, unsigned int D> 
    bool subtract(size_t N, const T* x, const T* y, T* r);

    /**
    * @brief multiply two vectors of values, r = x * y
    support in-place computation, e.g. x==r or y==r
    */
    template <typename T, unsigned int D> 
    bool multiply(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r);

    template <typename T, unsigned int D> 
    bool multiply(size_t N, const T* x, const T* y, T* r);

    /**
    * @brief divide two vectors of values, r = x / y
    support in-place computation, e.g. x==r
    no check for y==0
    */
    template <typename T, unsigned int D> 
    bool divide(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r);

    /**
    * @brief r = sqrt(x)
    */
    template <typename T, unsigned int D> 
    bool sqrt(const hoNDImage<T, D>& x, hoNDImage<T, D>& r);

    /**
    * @brief ind = min(abs(x(:))
    find the minimal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T, unsigned int D> 
    bool minAbsolute(const hoNDImage<T, D>& x, T& r, size_t& ind);

    /**
    * @brief ind = max(abs(x(:))
    find the miximal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T, unsigned int D> 
    bool maxAbsolute(const hoNDImage<T, D>& x, T& r, size_t& ind);

    /**
    * @brief r = x * conj(y)
    */
    template <typename T, unsigned int D> 
    bool multiplyConj(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r);

    /**
    * @brief r = conj(x)
    */
    template <typename T, unsigned int D> 
    bool conjugate(const hoNDImage<T, D>& x, hoNDImage<T, D>& r);

    /**
    * @brief if abs(x) is smaller than epsilon for its numeric type
    add epsilon to this x
    */
    template <typename T, unsigned int D> 
    bool addEpsilon(hoNDImage<T, D>& x);

    /**
    * @brief r = norm(x(:), 2)
    compute L2 norm of x
    */
    template <typename T, unsigned int D> 
    bool norm2(const hoNDImage<T, D>& x, typename realType<T>::Type& r);

    /**
    * @brief r = norm(x(:), 1)
    compute L1 norm of x = sum( abs(x(:) )
    */
    template <typename T, unsigned int D> 
    bool norm1(const hoNDImage<T, D>& x, typename realType<T>::Type& r);

    /**
    * @brief dot product of conj(x) and y
    r = conj(x) dot y
    */
    template <typename T, unsigned int D> 
    bool dotc(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, T& r);

    /**
    * @brief r = abs(x)
    */
    template <typename T, unsigned int D> 
    bool absolute(const hoNDImage<T, D>& x, hoNDImage<typename realType<T>::Type, D>& r);

    template <typename T, unsigned int D> 
    bool absolute(const hoNDImage< std::complex<T>, D >& x, hoNDImage< std::complex<T>, D >& r);

    /**
    * @brief r = angle(x)
    */
    template <typename T, unsigned int D> 
    bool argument(const hoNDImage<T, D>& x, hoNDImage<typename realType<T>::Type, D>& r);

    template <typename T, unsigned int D>
    bool argument(const hoNDImage<T, D>& x, hoNDImage< std::complex<T>, D>& r);

    /**
    * @brief r = 1/x
    */
    template <typename T, unsigned int D> 
    bool inv(const hoNDImage<T, D>& x, hoNDImage<T, D>& r);

    /**
    * @brief r = mean(x)
    */
    template <typename T, unsigned int D> 
    bool mean(const hoNDImage<T, D>& x, T& m);

    /**
    * @brief r = correlation_coefficient(a, b)
    */
    template <typename T, unsigned int D> 
    bool corrCoef(const hoNDImage<T, D>& a, const hoNDImage<T, D>& b, T& r);

    template<typename T, unsigned int D> void fill( hoNDImage<T, D>* x, T val );
    template<typename T, unsigned int D> void fill( hoNDImage<T, D>& x, T val );

    template<typename T, unsigned int D> void clear( hoNDImage<T, D>* x );
    template<typename T, unsigned int D> void clear( hoNDImage<T, D>& x );

    /**
    * @brief compute x *= a for an image
    */
    template <typename T, unsigned int D> bool scal(T a, hoNDImage<T, D>& x);
    template <typename T, unsigned int D> bool scal(T a, hoNDImage< std::complex<T>, D>& x);
}

#include "hoNDImage_util.hxx"
