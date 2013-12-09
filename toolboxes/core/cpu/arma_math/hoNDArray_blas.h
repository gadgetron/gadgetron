/** \file hoNDArray_blas.h
    \brief BLAS level-1 functions on the hoNDArray class.
    
    hoNDArray_blas.h provides BLAS level-1 functions on the hoNDArray class.
    The hoNDArray is temporarily reshaped to a column vector for the respective operations.
    The implementation is based on Armadillo.
    This code is purposely split into a header and underlying implementation (.cpp) 
    as this allows specific instantiation of the supported template types.     
    The supported types are float, double, std::complex<float>, std::complex<double>, 
    Gadgetron::complext<float>, and Gadgetron::complext<double>.
    There are currently no amin and amax functions instantiated for complex types 
    since Armadillo lacks an obvious method to compute the element-wise l1-norm.
*/

#pragma once

#include "hoNDArray.h"
#include "hoArmadillo.h"
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

namespace Gadgetron{

  /**
   * @brief Calculates the dot product of two arrays (as vectors).
   * @param[in] x Array 1. For complex arrays the complex conjugate of x is used.
   * @param[in] y Array 2.
   * @param[in] cc Specifies whether to use the complex conjugate of x (when applicable).
   * @return The dot product of x and y
   */
  template<class T> EXPORTCPUCOREMATH T dot( hoNDArray<T> *x, hoNDArray<T> *y, bool cc = true );

  /**
   * @brief Calculates the sum of the l1-norms of the array entries
   * @param[in] arr Input array
   * @return The l1-norm of the array
   */
  template<class T> EXPORTCPUCOREMATH typename realType<T>::Type asum( hoNDArray<T> *x );

  /**
   * @brief Calculates the sum of the l1-norms of the array entries
   * @param[in] arr Input array
   * @return The l1-norm of the array
   */
  template<class T> EXPORTCPUCOREMATH T asum( hoNDArray< std::complex<T> > *x );

  /**
   * @brief Calculates the sum of the l1-norms of the array entries
   * @param[in] arr Input array
   * @return The l1-norm of the array
   */
  template<class T> EXPORTCPUCOREMATH T asum( hoNDArray< complext<T> > *x );

  /**
   * @brief Calculates the l2-norm of the array (as a vector)
   * @param[in] arr Input array
   * @return The l2-norm of the array
   */
  template<class T> EXPORTCPUCOREMATH typename realType<T>::Type nrm2( hoNDArray<T> *x );

  /**
   * @brief Calculates the l1-norm of the array (as a vector)
   * @param[in] arr Input array
   * @return The l1-norm of the array
   */
  template<class T> EXPORTCPUCOREMATH typename realType<T>::Type nrm1( hoNDArray<T> *x );

  /**
   * @brief Returns the index of the array element with the smallest absolute value (l1 norm)
   * @param[in] x Input data
   * @return The array index corresponding to the smallest element in the array (0-indexing)
   */
  template<class T> EXPORTCPUCOREMATH unsigned long long amin( hoNDArray<T> *x );
 
  /**
   * @brief Returns the index of the array element with the smallest absolute value (l1 norm)
   * @param[in] x Input data
   * @return The array index corresponding to the smallest element in the array (0-indexing)
   */
  template<class T> EXPORTCPUCOREMATH unsigned long long amin( hoNDArray< std::complex<T> > *x );

  /**
   * @brief Returns the index of the array element with the smallest absolute value (l1 norm)
   * @param[in] x Input data
   * @return The array index corresponding to the smallest element in the array (0-indexing)
   */
  template<class T> EXPORTCPUCOREMATH unsigned long long amin( hoNDArray< complext<T> > *x );

  /**
   * @brief Returns the index of the array element with the largest absolute value (l1-norm)
   * @param[in] x Input data
   * @return The array index corresponding to the largest element in the array (0-indexing)
   */
  template<class T> EXPORTCPUCOREMATH unsigned long long amax( hoNDArray<T> *x );

  /**
   * @brief Returns the index of the array element with the largest absolute value (l1-norm)
   * @param[in] x Input data
   * @return The array index corresponding to the largest element in the array (0-indexing)
   */
  template<class T> EXPORTCPUCOREMATH unsigned long long amax( hoNDArray< std::complex<T> > *x );

  /**
   * @brief Returns the index of the array element with the largest absolute value (l1-norm)
   * @param[in] x Input data
   * @return The array index corresponding to the largest element in the array (0-indexing)
   */
  template<class T> EXPORTCPUCOREMATH unsigned long long amax( hoNDArray< complext<T> > *x );

  /**
   * @brief Calculates y = a*x+y in which x and y are considered as vectors
   * @param[in] a Scalar value
   * @param[in] x Array
   * @param[in,out] y Array
   */
  template<class T> EXPORTCPUCOREMATH void axpy( T a, hoNDArray<T> *x, hoNDArray<T> *y );

    /**
     * Besides the functions calling the arma, there are some more functions directly calling the MKL routines
     */

    #ifdef USE_MKL

    template<> EXPORTCPUCOREMATH float nrm1( hoNDArray<float> *x );
    template<> EXPORTCPUCOREMATH double nrm1( hoNDArray<double> *x );

    // BLAS dotc and dotu
    // res = conj(x) dot y
    EXPORTCPUCOREMATH GT_Complex8 dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y);
    EXPORTCPUCOREMATH GT_Complex16 dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y);

    // res = x dot y
    EXPORTCPUCOREMATH GT_Complex8 dotu(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y);
    EXPORTCPUCOREMATH GT_Complex16 dotu(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y);

    // other variants for axpy
    // r = a*x+y
    EXPORTCPUCOREMATH bool axpy(float a, const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    EXPORTCPUCOREMATH bool axpy(double a, const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool axpy(const GT_Complex8& a, const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
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

    // sort the vector
    // isascending: true for ascending and false for descending
    EXPORTCPUCOREMATH bool sort(const hoNDArray<float>& x, hoNDArray<float>& r, bool isascending);
    EXPORTCPUCOREMATH bool sort(const hoNDArray<double>& x, hoNDArray<double>& r, bool isascending);

    #endif // USE_MKL
}
