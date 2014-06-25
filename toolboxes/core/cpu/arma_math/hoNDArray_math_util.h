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

#include "complext.h"
#include "cpucore_math_export.h"
#include "GadgetronCommon.h"
#include <complex>

#ifdef USE_MKL
#include "mkl.h"
#endif // USE_MKL

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
#ifndef USE_MKL

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
    * @brief compute x *= a for an image
    */
    template <typename T, unsigned int D> bool scal(T a, hoNDImage<T, D>& x);
    template <typename T, unsigned int D> bool scal(T a, hoNDImage< std::complex<T>, D>& x);

    /**
    * @brief sort the ND array
    */
    template <typename T> bool sort(const hoNDArray<T>& x, hoNDArray<T>& r, bool isascending);

#endif // USE_MKL

#ifdef USE_MKL

    // besides the arma calls, some functions are implemented with the MKL vector utilities

    // ------------------------------------------------------------------
    // float
    // ------------------------------------------------------------------
    EXPORTCPUCOREMATH bool add(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r); // r = x + y
    EXPORTCPUCOREMATH bool add(size_t N, const float* x, const float* y, float* r);
    EXPORTCPUCOREMATH bool subtract(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r); // r = x - y
    EXPORTCPUCOREMATH bool multiply(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r); // r = x * y
    EXPORTCPUCOREMATH bool divide(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r); // r = x / y
    EXPORTCPUCOREMATH bool absolute(const hoNDArray<float>& x, hoNDArray<float>& r); // r = abs(x)
    EXPORTCPUCOREMATH bool argument(const hoNDArray<float>& x, hoNDArray<float>& r); // r = angle(x)
    EXPORTCPUCOREMATH bool sqrt(const hoNDArray<float>& x, hoNDArray<float>& r); // r = sqrt(x)
    EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<float>& x, float& r, size_t& ind); // minimal absolute value and index
    EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<float>& x, float& r, size_t& ind); // maximal absolute value and index
    EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<float>& x); // x = x + Epsilon if x==0, prepare for division
    EXPORTCPUCOREMATH bool norm2(const hoNDArray<float>& x, float& r);
    EXPORTCPUCOREMATH bool norm1(const hoNDArray<float>& x, float& r);
    EXPORTCPUCOREMATH bool conv2(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& z); // x: input data, y: convolution kernel, z: output; each 2D slice is convolved
    EXPORTCPUCOREMATH bool conv3(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& z); // x: input data, y: convolution kernel, z: output; each 3D volume is convolved
    EXPORTCPUCOREMATH bool inv(const hoNDArray<float>& x, hoNDArray<float>& r); // r = 1/x
    EXPORTCPUCOREMATH bool axpy(float a, const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r); // r = a*x+y

    // ------------------------------------------------------------------
    // double
    // ------------------------------------------------------------------
    EXPORTCPUCOREMATH bool add(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool add(size_t N, const double* x, const double* y, double* r);
    EXPORTCPUCOREMATH bool subtract(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool multiply(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool divide(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool absolute(const hoNDArray<double>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool argument(const hoNDArray<double>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool sqrt(const hoNDArray<double>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<double>& x, double& r, size_t& ind);
    EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<double>& x, double& r, size_t& ind);
    EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<double>& x);
    EXPORTCPUCOREMATH bool norm2(const hoNDArray<double>& x, double& r);
    EXPORTCPUCOREMATH bool norm1(const hoNDArray<double>& x, double& r);
    EXPORTCPUCOREMATH bool conv2(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& z);
    EXPORTCPUCOREMATH bool conv3(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& z);
    EXPORTCPUCOREMATH bool inv(const hoNDArray<double>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool axpy(double a, const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);

    // ------------------------------------------------------------------
    // complex float
    // ------------------------------------------------------------------
    EXPORTCPUCOREMATH bool add(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool add(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r);
    EXPORTCPUCOREMATH bool subtract(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool subtract(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r);
    EXPORTCPUCOREMATH bool multiply(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool multiply(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r);
    EXPORTCPUCOREMATH bool divide(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r);
    EXPORTCPUCOREMATH bool sqrt(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind);
    EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind);
    EXPORTCPUCOREMATH bool multiplyConj(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool argument(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r);
    EXPORTCPUCOREMATH bool conjugate(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<GT_Complex8>& x);
    EXPORTCPUCOREMATH bool norm2(const hoNDArray<GT_Complex8>& x, float& r);
    EXPORTCPUCOREMATH bool norm1(const hoNDArray<GT_Complex8>& x, float& r);
    EXPORTCPUCOREMATH bool dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, GT_Complex8& r);
    EXPORTCPUCOREMATH GT_Complex8 dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y);
    EXPORTCPUCOREMATH GT_Complex8 dotu(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y); // res = x dot y

    EXPORTCPUCOREMATH bool conv2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& z);
    EXPORTCPUCOREMATH bool conv3(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& z);

    EXPORTCPUCOREMATH bool corr2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& z); // x: input data [RO E1 ...], y: corr kernel [kro ke1], z: output; each 2D slice is correlated
    EXPORTCPUCOREMATH bool corr3(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& z); // x: input data [RO E1 E2 ...], y: corr kernel [kro ke1 ke2], z: output; each 3D volume is correlated

    EXPORTCPUCOREMATH bool inv(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool axpy(const GT_Complex8& a, const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);

    // ------------------------------------------------------------------
    // complex double
    // ------------------------------------------------------------------
    EXPORTCPUCOREMATH bool add(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool add(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r);
    EXPORTCPUCOREMATH bool subtract(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool subtract(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r);
    EXPORTCPUCOREMATH bool multiply(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool multiply(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r);
    EXPORTCPUCOREMATH bool divide(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool sqrt(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind);
    EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind);
    EXPORTCPUCOREMATH bool multiplyConj(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool argument(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool conjugate(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<GT_Complex16>& x);
    EXPORTCPUCOREMATH bool norm2(const hoNDArray<GT_Complex16>& x, double& r);
    EXPORTCPUCOREMATH bool norm1(const hoNDArray<GT_Complex16>& x, double& r);
    EXPORTCPUCOREMATH bool dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, GT_Complex16& r);
    EXPORTCPUCOREMATH GT_Complex16 dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y);
    EXPORTCPUCOREMATH GT_Complex16 dotu(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y);

    EXPORTCPUCOREMATH bool conv2(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& z);
    EXPORTCPUCOREMATH bool conv3(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& z);
    EXPORTCPUCOREMATH bool corr2(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& z);
    EXPORTCPUCOREMATH bool corr3(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& z);

    EXPORTCPUCOREMATH bool inv(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool axpy(const GT_Complex16& a, const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    // vector-scalar product
    // r = a*x
    EXPORTCPUCOREMATH bool scal(float a, hoNDArray<float>& x);
    EXPORTCPUCOREMATH bool scal(double a, hoNDArray<double>& x);
    EXPORTCPUCOREMATH bool scal(float a, hoNDArray<GT_Complex8>& x);
    EXPORTCPUCOREMATH bool scal(double a, hoNDArray<GT_Complex16>& x);
    EXPORTCPUCOREMATH bool scal(GT_Complex8 a, hoNDArray<GT_Complex8>& x);
    EXPORTCPUCOREMATH bool scal(GT_Complex16 a, hoNDArray<GT_Complex16>& x);

    EXPORTCPUCOREMATH bool scal(float a, float*x, long long N);
    EXPORTCPUCOREMATH bool scal(double a, double*x, long long N);
    EXPORTCPUCOREMATH bool scal(float a, GT_Complex8*x, long long N);
    EXPORTCPUCOREMATH bool scal(double a, GT_Complex16*x, long long N);
    EXPORTCPUCOREMATH bool scal(GT_Complex8 a, GT_Complex8*x, long long N);
    EXPORTCPUCOREMATH bool scal(GT_Complex16 a, GT_Complex16*x, long long N);

    template <unsigned int D> EXPORTCPUCOREMATH bool scal(float a, hoNDImage<float, D>& x);
    template <unsigned int D> EXPORTCPUCOREMATH bool scal(double a, hoNDImage<double, D>& x);
    template <unsigned int D> EXPORTCPUCOREMATH bool scal(float a, hoNDImage<GT_Complex8, D>& x);
    template <unsigned int D> EXPORTCPUCOREMATH bool scal(double a, hoNDImage<GT_Complex16, D>& x);
    template <unsigned int D> EXPORTCPUCOREMATH bool scal(GT_Complex8 a, hoNDImage<GT_Complex8, D>& x);
    template <unsigned int D> EXPORTCPUCOREMATH bool scal(GT_Complex16 a, hoNDImage<GT_Complex16, D>& x);

    // sort the vector
    // isascending: true for ascending and false for descending
    EXPORTCPUCOREMATH bool sort(const hoNDArray<float>& x, hoNDArray<float>& r, bool isascending);
    EXPORTCPUCOREMATH bool sort(const hoNDArray<double>& x, hoNDArray<double>& r, bool isascending);

#endif // USE_MKL

    // ------------------------------------------------------------------

    // fill in an array
    template<typename T> void fill( size_t N, T* pX, T val );

    template<typename T, unsigned int D> void fill( hoNDImage<T, D>* x, T val );
    template<typename T, unsigned int D> void fill( hoNDImage<T, D>& x, T val );

    template<typename T, unsigned int D> void clear( hoNDImage<T, D>* x );
    template<typename T, unsigned int D> void clear( hoNDImage<T, D>& x );

    // r = abs(x)
    template <typename T> 
    bool absolute(const hoNDArray< std::complex<T> >& x, hoNDArray< std::complex<T> >& r);
}

#include "hoNDArray_math_util.hxx"
