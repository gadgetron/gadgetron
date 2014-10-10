/** \file  hoNDArray_math_util.h
    \brief math functions for hoNDArray and hoNDImage not using armadillo
           Some functions are using MKL; some have the general implementation
           The MKL implementation has priority because of its speed
*/

#pragma once

#include "hoNDArray.h"
#include "ho2DArray.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include "ho5DArray.h"
#include "ho6DArray.h"
#include "ho7DArray.h"
#include "hoNDImage.h"
#include "hoNDMath_util.h"

#include "complext.h"
#include "cpucore_math_export.h"
#include "GadgetronCommon.h"
#include <complex>

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
    * @brief add two vectors of values, r = x + y
    support in-place computation, e.g. x==r or y==r
    */
    template <typename T> 
    bool add(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief subtract two vectors of values, r = x - y
    support in-place computation, e.g. x==r
    */
    template <typename T> 
    bool subtract(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    template <typename T> 
    bool subtract(size_t N, const T* x, const T* y, T* r);

    /**
    * @brief multiply two vectors of values, r = x * y
    support in-place computation, e.g. x==r or y==r
    */
    template <typename T> 
    bool multiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    template <typename T> 
    bool multiply(size_t N, const T* x, const T* y, T* r);

    /**
    * @brief divide two vectors of values, r = x / y
    support in-place computation, e.g. x==r
    no check for y==0
    */
    template <typename T> 
    bool divide(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief r = sqrt(x)
    */
    template <typename T> 
    bool sqrt(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief ind = min(abs(x(:))
    find the minimal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T> 
    bool minAbsolute(const hoNDArray<T>& x, T& r, size_t& ind);

    /**
    * @brief ind = max(abs(x(:))
    find the miximal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T> 
    bool maxAbsolute(const hoNDArray<T>& x, T& r, size_t& ind);

    /**
    * @brief r = x * conj(y)
    */
    template <typename T> 
    bool multiplyConj(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief r = conj(x)
    */
    template <typename T> 
    bool conjugate(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief if abs(x) is smaller than epsilon for its numeric type
    add epsilon to this x
    */
    template <typename T> 
    bool addEpsilon(hoNDArray<T>& x);

    /**
    * @brief r = norm(x(:), 2)
    compute L2 norm of x
    */
    template <typename T> 
    bool norm2(const hoNDArray<T>& x, typename realType<T>::Type& r);

    /**
    * @brief r = norm(x(:), 1)
    compute L1 norm of x = sum( abs(x(:) )
    */
    template <typename T> 
    bool norm1(const hoNDArray<T>& x, typename realType<T>::Type& r);

    /**
    * @brief dot product of conj(x) and y
    r = conj(x) dot y
    */
    template <typename T> 
    bool dotc(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r);

    /**
    * @brief r = abs(x)
    */
    template <typename T> 
    bool absolute(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r);

    /**
    * @brief r = angle(x)
    */
    template <typename T> 
    bool argument(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r);

    /**
    * @brief r = 1/x
    */
    template <typename T> 
    bool inv(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief 2D convolution
             x: input data, y: convolution kernel, z: output; each 2D slice is convolved
    */
    template<typename T> 
    bool conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z);

    /**
    * @brief 3D convolution
             x: input data, y: convolution kernel, z: output; each 3D volume is convolved
    */
    template<typename T> 
    bool conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z);

    /**
    * @brief compute conj(x) dot y
    */
    template <typename T> T dotc(const hoNDArray<T>& x, const hoNDArray<T>& y);

    /**
    * @brief compute x dot y
    */
    template <typename T> T dotu(const hoNDArray<T>& x, const hoNDArray<T>& y);

    /**
    * @brief compute r = a*x + y
    */
    template <typename T> bool axpy(T a, const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief compute x *= a
    */
    template <typename T> bool scal(T a, hoNDArray<T>& x);
    template <typename T> bool scal(T a, hoNDArray< std::complex<T> >& x);

    /**
    * @brief compute x *= a for an array
    */
    template <typename T> bool scal(T a, T*x, long long N);
    template <typename T> bool scal(T a, std::complex<T>*x, long long N);

    /**
    * @brief sort the ND array
    */
    template <typename T> bool sort(const hoNDArray<T>& x, hoNDArray<T>& r, bool isascending);

    // ------------------------------------------------------------------

    // fill in an array
    template<typename T> void fill( size_t N, T* pX, T val );

    // r = abs(x)
    template <typename T> 
    bool absolute(const hoNDArray< std::complex<T> >& x, hoNDArray< std::complex<T> >& r);
}

#include "hoNDArray_math_util.hxx"
