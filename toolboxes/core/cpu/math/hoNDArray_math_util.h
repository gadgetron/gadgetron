/** \file  hoNDArray_math_util.h
    \brief math functions for hoNDArray and hoNDImage not using armadillo
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

#include "complext.h"
#include "cpucore_math_export.h"
#include "GadgetronCommon.h"
#include <complex>

namespace Gadgetron
{

    /**
    * @brief add two vectors of values, r = x + y
    support in-place computation, e.g. x==r or y==r
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool add(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief subtract two vectors of values, r = x - y
    support in-place computation, e.g. x==r
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool subtract(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief multiply two vectors of values, r = x * y
    support in-place computation, e.g. x==r or y==r
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool multiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief divide two vectors of values, r = x / y
    support in-place computation, e.g. x==r
    no check for y==0
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool divide(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief r = sqrt(x)
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool sqrt(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief ind = min(abs(x(:))
    find the minimal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool minAbsolute(const hoNDArray<T>& x, T& r, size_t& ind);

    /**
    * @brief ind = max(abs(x(:))
    find the miximal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool maxAbsolute(const hoNDArray<T>& x, T& r, size_t& ind);

    /**
    * @brief r = x * conj(y)
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool multiplyConj(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief r = conj(x)
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool conjugate(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief if abs(x) is smaller than epsilon for its numeric type
    add epsilon to this x
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool addEpsilon(hoNDArray<T>& x);

    /**
    * @brief r = norm(x(:), 2)
    compute L2 norm of x
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool norm2(const hoNDArray<T>& x, typename realType<T>::Type& r);

    template <typename T> EXPORTCPUCOREMATH 
    typename realType<T>::Type norm2(const hoNDArray<T>& x);

    /**
    * @brief r = norm(x(:), 1)
    compute L1 norm of x = sum( abs(x(:) )
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool norm1(const hoNDArray<T>& x, typename realType<T>::Type& r);

    template <typename T> EXPORTCPUCOREMATH 
    typename realType<T>::Type norm1(const hoNDArray<T>& x);

    /**
    * @brief dot product of conj(x) and y
    r = conj(x) dot y
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool dotc(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r);

    template <typename T> EXPORTCPUCOREMATH 
    T dotc(const hoNDArray<T>& x, const hoNDArray<T>& y);

    /**
    * @brief dot product of x and y
    r = x dot y
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool dotu(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r);

    template <typename T> EXPORTCPUCOREMATH 
    T dotu(const hoNDArray<T>& x, const hoNDArray<T>& y);

    /**
    * @brief r = abs(x)
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool absolute(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r);

    /**
    * @brief absolute of a complex array
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool absolute(const hoNDArray< std::complex<T> >& x, hoNDArray< std::complex<T> >& r);

    /**
    * @brief r = angle(x)
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool argument(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r);

    /**
    * @brief r = 1/x
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool inv(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief 2D convolution
             x: input data, y: convolution kernel, z: output; each 2D slice is convolved
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z);

    /**
    * @brief 3D convolution
             x: input data, y: convolution kernel, z: output; each 3D volume is convolved
    */
    template <typename T> EXPORTCPUCOREMATH 
    bool conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z);

    /**
    * @brief compute r = a*x + y
    */
    template <typename T> EXPORTCPUCOREMATH bool axpy(T a, const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief compute x *= a
    */
    template <typename T> EXPORTCPUCOREMATH bool scal(T a, hoNDArray<T>& x);
    template <typename T> EXPORTCPUCOREMATH bool scal(T a, hoNDArray< std::complex<T> >& x);

    /**
    * @brief compute x *= a for an array
    */
    template <typename T> EXPORTCPUCOREMATH bool scal(T a, T*x, long long N);
    template <typename T> EXPORTCPUCOREMATH bool scal(T a, std::complex<T>*x, long long N);

    /**
    * @brief sort the ND array
    */
    template <typename T> EXPORTCPUCOREMATH bool sort(const hoNDArray<T>& x, hoNDArray<T>& r, bool isascending);

    /**
    * @brief fill in an array
    */
    template <typename T> EXPORTCPUCOREMATH void fill( hoNDArray<T>* x, T val);
    template <typename T> EXPORTCPUCOREMATH void fill( hoNDArray<T>& x, T val );

    /**
    * @brief computes the sum of magnitudes of the vector elements.
    */
    template<class T> EXPORTCPUCOREMATH void asum(const hoNDArray<T>& x, typename realType<T>::Type& r);
    template<class T> EXPORTCPUCOREMATH typename realType<T>::Type asum(const hoNDArray<T>& x);

    /**
    * @brief finds the index of the element with the maximal absolute value.
    */
    template<class T> EXPORTCPUCOREMATH size_t amax(const hoNDArray<T>& x);

    /**
    * @brief real and imag to complex
    */
    template<class T> EXPORTCPUCOREMATH 
    bool real_imag_to_complex(const hoNDArray<typename realType<T>::Type>& real, const hoNDArray<typename realType<T>::Type>& imag, hoNDArray<T>& cplx);

    /**
    * @brief complex to real and imag
    */
    template<class T> EXPORTCPUCOREMATH 
    bool complex_to_real_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real, hoNDArray<typename realType<T>::Type>& imag);

    /**
    * @brief get the real part of complex
    */
    template<class T> EXPORTCPUCOREMATH 
    bool complex_to_real(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real);

    template<class T> EXPORTCPUCOREMATH 
    bool complex_to_real(const hoNDArray<T>& cplx, hoNDArray<T>& real);

    template<class T> 
    bool complex_to_real(hoNDArray<T>& cplx);

    /**
    * @brief get the imag part of complex
    */
    template<class T> EXPORTCPUCOREMATH 
    bool complex_to_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& imag);

    template<class T> EXPORTCPUCOREMATH 
    bool complex_to_imag(const hoNDArray<T>& cplx, hoNDArray<T>& imag);

    template<class T> EXPORTCPUCOREMATH 
    bool complex_to_imag(hoNDArray<T>& cplx);

    /**
    * @brief get complex array whose real part is the input and imag part is zero
    */
    template<class T> EXPORTCPUCOREMATH 
    bool real_to_complex(const hoNDArray<typename realType<T>::Type>& real, hoNDArray<T>& cplx);

    /**
    * @brief get the min and max value from an array (only for float and double type)
    */
    template <class T> EXPORTCPUCOREMATH 
    bool minValue(const hoNDArray<T>& a, T& v);

    template <class T> EXPORTCPUCOREMATH 
    bool maxValue(const hoNDArray<T>& a, T& v);
}
