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

#ifdef USE_MKL
#include "mkl.h"
#endif // USE_MKL

#include "hoNDArray_math_util.h"
#include "hoNDInterpolator.h"

#ifdef GT_Complex8
#undef GT_Complex8
#endif // GT_Complex8
typedef std::complex<float> GT_Complex8;

#ifdef GT_Complex16
#undef GT_Complex16
#endif // GT_Complex16
typedef std::complex<double> GT_Complex16;

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

    /// compute the gradient for an ND image
    /// the central difference is computed, the border-value boundary condition is used
    template<class T, unsigned int D> EXPORTCPUCOREMATH bool gradient(const hoNDImage<T, D>& x, hoNDImage<T, D> gx[]);

    /// perform the gaussian filter for every dimension
    /// sigma is in the unit of pixel
    template<class ArrayType, class T2> EXPORTCPUCOREMATH bool filterGaussian(ArrayType& x, T2 sigma[]);

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

#ifdef USE_MKL

    // r = x + y
    template <unsigned int D> EXPORTCPUCOREMATH bool add(const hoNDImage<float, D>& x, const hoNDImage<float, D>& y, hoNDImage<float, D>& r);
    // r = x - y
    template <unsigned int D> EXPORTCPUCOREMATH bool subtract(const hoNDImage<float, D>& x, const hoNDImage<float, D>& y, hoNDImage<float, D>& r);
    // r = x * y
    template <unsigned int D> EXPORTCPUCOREMATH bool multiply(const hoNDImage<float, D>& x, const hoNDImage<float, D>& y, hoNDImage<float, D>& r);
    // r = x / y
    template <unsigned int D> EXPORTCPUCOREMATH bool divide(const hoNDImage<float, D>& x, const hoNDImage<float, D>& y, hoNDImage<float, D>& r);
    // r = abs(x)
    template <unsigned int D> EXPORTCPUCOREMATH bool absolute(const hoNDImage<float, D>& x, hoNDImage<float, D>& r);
    // r = angle(x)
    template <unsigned int D> EXPORTCPUCOREMATH bool argument(const hoNDImage<float, D>& x, hoNDImage<float, D>& r);
    // r = sqrt(x)
    template <unsigned int D> EXPORTCPUCOREMATH bool sqrt(const hoNDImage<float, D>& x, hoNDImage<float, D>& r);
    // minimal absolute value and index
    template <unsigned int D> EXPORTCPUCOREMATH bool minAbsolute(const hoNDImage<float, D>& x, float& r, size_t& ind);
    // maximal absolute value and index
    template <unsigned int D> EXPORTCPUCOREMATH bool maxAbsolute(const hoNDImage<float, D>& x, float& r, size_t& ind);
    // x = x + Epsilon if x==0, prepare for division
    template <unsigned int D> EXPORTCPUCOREMATH bool addEpsilon(hoNDImage<float, D>& x);
    template <unsigned int D> EXPORTCPUCOREMATH bool norm2(const hoNDImage<float, D>& x, float& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool norm1(const hoNDImage<float, D>& x, float& r);
    // x: input data, y: convolution kernel, z: output; each 2D slice is convolved
    template <unsigned int D> EXPORTCPUCOREMATH bool conv2(const hoNDImage<float, D>& x, const hoNDImage<float, D>& y, hoNDImage<float, D>& z);
    // x: input data, y: convolution kernel, z: output; each 3D volume is convolved
    template <unsigned int D> EXPORTCPUCOREMATH bool conv3(const hoNDImage<float, D>& x, const hoNDImage<float, D>& y, hoNDImage<float, D>& z);
    // r = 1/x
    template <unsigned int D> EXPORTCPUCOREMATH bool inv(const hoNDImage<float, D>& x, hoNDImage<float, D>& r);

    template <unsigned int D> EXPORTCPUCOREMATH bool add(const hoNDImage<double, D>& x, const hoNDImage<double, D>& y, hoNDImage<double, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool subtract(const hoNDImage<double, D>& x, const hoNDImage<double, D>& y, hoNDImage<double, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool multiply(const hoNDImage<double, D>& x, const hoNDImage<double, D>& y, hoNDImage<double, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool divide(const hoNDImage<double, D>& x, const hoNDImage<double, D>& y, hoNDImage<double, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool absolute(const hoNDImage<double, D>& x, hoNDImage<double, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool argument(const hoNDImage<double, D>& x, hoNDImage<double, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool sqrt(const hoNDImage<double, D>& x, hoNDImage<double, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool minAbsolute(const hoNDImage<double, D>& x, double& r, size_t& ind);
    template <unsigned int D> EXPORTCPUCOREMATH bool maxAbsolute(const hoNDImage<double, D>& x, double& r, size_t& ind);
    template <unsigned int D> EXPORTCPUCOREMATH bool addEpsilon(hoNDImage<double, D>& x);
    template <unsigned int D> EXPORTCPUCOREMATH bool norm2(const hoNDImage<double, D>& x, double& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool norm1(const hoNDImage<double, D>& x, double& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool conv2(const hoNDImage<double, D>& x, const hoNDImage<double, D>& y, hoNDImage<double, D>& z);
    template <unsigned int D> EXPORTCPUCOREMATH bool conv3(const hoNDImage<double, D>& x, const hoNDImage<double, D>& y, hoNDImage<double, D>& z);
    template <unsigned int D> EXPORTCPUCOREMATH bool inv(const hoNDImage<double, D>& x, hoNDImage<double, D>& r);

    template <unsigned int D> EXPORTCPUCOREMATH bool add(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool subtract(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool multiply(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool divide(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool absolute(const hoNDImage<GT_Complex8, D>& x, hoNDImage<float, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool absolute(const hoNDImage<GT_Complex8, D>& x, hoNDImage<GT_Complex8, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool sqrt(const hoNDImage<GT_Complex8, D>& x, hoNDImage<GT_Complex8, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool minAbsolute(const hoNDImage<GT_Complex8, D>& x, GT_Complex8& r, size_t& ind);
    template <unsigned int D> EXPORTCPUCOREMATH bool maxAbsolute(const hoNDImage<GT_Complex8, D>& x, GT_Complex8& r, size_t& ind);
    template <unsigned int D> EXPORTCPUCOREMATH bool multiplyConj(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& r); // r = x * conj(y)
    template <unsigned int D> EXPORTCPUCOREMATH bool argument(const hoNDImage<GT_Complex8, D>& x, hoNDImage<float, D>& r); // r = angle(x)
    template <unsigned int D> EXPORTCPUCOREMATH bool conjugate(const hoNDImage<GT_Complex8, D>& x, hoNDImage<GT_Complex8, D>& r); // r = conj(x)
    template <unsigned int D> EXPORTCPUCOREMATH bool addEpsilon(hoNDImage<GT_Complex8, D>& x);
    template <unsigned int D> EXPORTCPUCOREMATH bool norm2(const hoNDImage<GT_Complex8, D>& x, float& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool norm1(const hoNDImage<GT_Complex8, D>& x, float& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool dotc(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, GT_Complex8& r); // x'*y, x and y are N*1 vector
    template <unsigned int D> EXPORTCPUCOREMATH bool conv2(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& z);
    template <unsigned int D> EXPORTCPUCOREMATH bool conv3(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& z);
    template <unsigned int D> EXPORTCPUCOREMATH bool corr2(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& z); // x: input data [RO E1 ...], y: corr kernel [kro ke1], z: output; each 2D slice is correlated
    template <unsigned int D> EXPORTCPUCOREMATH bool corr3(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& z); // x: input data [RO E1 E2 ...], y: corr kernel [kro ke1 ke2], z: output; each 3D volume is correlated
    template <unsigned int D> EXPORTCPUCOREMATH bool inv(const hoNDImage<GT_Complex8, D>& x, hoNDImage<GT_Complex8, D>& r);

    template <unsigned int D> EXPORTCPUCOREMATH bool add(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool subtract(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool multiply(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool divide(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool absolute(const hoNDImage<GT_Complex16, D>& x, hoNDImage<double, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool absolute(const hoNDImage<GT_Complex16, D>& x, hoNDImage<GT_Complex16, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool sqrt(const hoNDImage<GT_Complex16, D>& x, hoNDImage<GT_Complex16, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool minAbsolute(const hoNDImage<GT_Complex16, D>& x, GT_Complex16& r, size_t& ind);
    template <unsigned int D> EXPORTCPUCOREMATH bool maxAbsolute(const hoNDImage<GT_Complex16, D>& x, GT_Complex16& r, size_t& ind);
    template <unsigned int D> EXPORTCPUCOREMATH bool multiplyConj(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool argument(const hoNDImage<GT_Complex16, D>& x, hoNDImage<double, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool conjugate(const hoNDImage<GT_Complex16, D>& x, hoNDImage<GT_Complex16, D>& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool addEpsilon(hoNDImage<GT_Complex16, D>& x);
    template <unsigned int D> EXPORTCPUCOREMATH bool norm2(const hoNDImage<GT_Complex16, D>& x, double& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool norm1(const hoNDImage<GT_Complex16, D>& x, double& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool dotc(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, GT_Complex16& r);
    template <unsigned int D> EXPORTCPUCOREMATH bool conv2(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& z);
    template <unsigned int D> EXPORTCPUCOREMATH bool conv3(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& z);
    template <unsigned int D> EXPORTCPUCOREMATH bool corr2(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& z);
    template <unsigned int D> EXPORTCPUCOREMATH bool corr3(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& z);
    template <unsigned int D> EXPORTCPUCOREMATH bool inv(const hoNDImage<GT_Complex16, D>& x, hoNDImage<GT_Complex16, D>& r);

#endif // USE_MKL
}

#include "hoNDImage_util.hxx"
