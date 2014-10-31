/** \file  hoNDMath_util.h
    \brief math function utility
*/

#pragma once

#include "cpucore_math_export.h"
#include <complex>
#include "hoNDArray.h"

namespace Gadgetron { namespace math { 

    /**
    * @brief compute r = a*x + y
    */
    template <typename T> EXPORTCPUCOREMATH void axpy(T a, size_t N, const T* x, const T* y, T* r);

    /**
    * @brief add two vectors of values, r = x + y
    support in-place computation, e.g. x==r or y==r
    */
    template <typename T> EXPORTCPUCOREMATH void add(size_t N, const T* x, const T* y, T* r);

    /**
    * @brief subtract two vectors of values, r = x - y
    support in-place computation, e.g. x==r
    */
    template <typename T> EXPORTCPUCOREMATH void subtract(size_t N, const T* x, const T* y, T* r);

    /**
    * @brief multiply two vectors of values, r = x * y
    support in-place computation, e.g. x==r or y==r
    */
    template <typename T> EXPORTCPUCOREMATH void multiply(size_t N, const T* x, const T* y, T* r);

    /**
    * @brief divide two vectors of values, r = x / y
    support in-place computation, e.g. x==r
    no check for y==0
    */
    template <typename T> EXPORTCPUCOREMATH void divide(size_t N, const T* x, const T* y, T* r);

    /**
    * @brief r = sqrt(x)
    */
    template <typename T> EXPORTCPUCOREMATH void sqrt(size_t N, const T* x, T* r);

    /**
    * @brief ind = min(abs(x(:))
    find the minimal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T> EXPORTCPUCOREMATH void minAbsolute(size_t N, const T* x, T& r, size_t& ind);

    /**
    * @brief ind = max(abs(x(:))
    find the miximal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T> EXPORTCPUCOREMATH void maxAbsolute(size_t N, const T* x, T& r, size_t& ind);

    /**
    * @brief r = x * conj(y)
    */
    template <typename T> EXPORTCPUCOREMATH void multiplyConj(size_t N, const T* x, const T* y, T* r);

    /**
    * @brief r = conj(x)
    */
    template <typename T> EXPORTCPUCOREMATH void conjugate(size_t N, const T* x, T* r);

    /**
    * @brief if abs(x) is smaller than epsilon for its numeric type
    add epsilon to this x
    */
    template <typename T> EXPORTCPUCOREMATH void addEpsilon(size_t N, T* x);

    /**
    * @brief r = norm(x(:), 2)
    compute L2 norm of x
    */
    template <typename T> EXPORTCPUCOREMATH void norm2(size_t N, const T* x, typename realType<T>::Type& r);
    template <typename T> EXPORTCPUCOREMATH typename realType<T>::Type norm2(size_t N, const T* x);

    /**
    * @brief r = norm(x(:), 1)
    compute L1 norm of x = sum( abs(x(:) )
    */
    template <typename T> EXPORTCPUCOREMATH void norm1(size_t N, const T* x, typename realType<T>::Type& r);
    template <typename T> EXPORTCPUCOREMATH typename realType<T>::Type norm1(size_t N, const T* x);

    /**
    * @brief dot product of conj(x) and y
    r = conj(x) dot y
    */
    template <typename T> EXPORTCPUCOREMATH void dotc(size_t N, const T* x, const T* y, T& r);
    template <typename T> EXPORTCPUCOREMATH T dotc(size_t N, const T* x, const T* y);

    /**
    * @brief compute x dot y
    */
    template <typename T> EXPORTCPUCOREMATH void dotu(size_t N, const T* x, const T* y, T& r);
    template <typename T> EXPORTCPUCOREMATH T dotu(size_t N, const T* x, const T* y);

    /**
    * @brief computes the sum of magnitudes of the vector elements.
    */
    template<class T> EXPORTCPUCOREMATH void asum(size_t N, const T* x, typename realType<T>::Type& r);
    template<class T> EXPORTCPUCOREMATH typename realType<T>::Type asum(size_t N, const T* x);

    /**
    * @brief finds the index of the element with the maximal absolute value.
    */
    template<class T> EXPORTCPUCOREMATH size_t amax(size_t N, const T* x);

    /**
    * @brief r = abs(x)
    */
    template <typename T> EXPORTCPUCOREMATH void absolute(size_t N, const T* x, typename realType<T>::Type* r);

    template <typename T> EXPORTCPUCOREMATH void absolute(size_t N, const std::complex<T>* x, std::complex<T>* r);

    /**
    * @brief r = argument(x), angle of x
    */
    template <typename T> EXPORTCPUCOREMATH void argument(size_t N, const T* x, typename realType<T>::Type* r);

    /**
    * @brief r = 1/x
    */
    template <typename T> EXPORTCPUCOREMATH void inv(size_t N, const T* x, T* r);

    /**
    * @brief 2D convolution
             x: input data, y: convolution kernel, z: output; each 2D slice is convolved
    */
    template<typename T> EXPORTCPUCOREMATH void conv2(size_t RO, size_t E1, size_t num, const T* x, size_t kRO, size_t kE1, const T* y, T* z);

    /**
    * @brief 3D convolution
             x: input data, y: convolution kernel, z: output; each 3D volume is convolved
    */
    template<typename T> EXPORTCPUCOREMATH void conv3(size_t RO, size_t E1, size_t E2, size_t num, const T* x, size_t kRO, size_t kE1, size_t kE2, const T* y, T* z);

    /**
    * @brief compute x *= a for an array
    */
    template <typename T> EXPORTCPUCOREMATH void scal(size_t N, T a, T*x);
    template <typename T> EXPORTCPUCOREMATH void scal(size_t N, T a, std::complex<T>*x);

    /**
    * @brief sort the ND array
    */
    template <typename T> EXPORTCPUCOREMATH void sort(size_t N, const T* x, T* r, bool isascending);
                    
    /**
    * @brief fill in an array
    */
    template<typename T> EXPORTCPUCOREMATH void fill( size_t N, T* x, T v );
}}
