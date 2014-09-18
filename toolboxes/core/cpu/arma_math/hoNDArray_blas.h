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
    template<class T> EXPORTCPUCOREMATH size_t amin( hoNDArray<T> *x );

    /**
    * @brief Returns the index of the array element with the smallest absolute value (l1 norm)
    * @param[in] x Input data
    * @return The array index corresponding to the smallest element in the array (0-indexing)
    */
    template<class T> EXPORTCPUCOREMATH size_t amin( hoNDArray< std::complex<T> > *x );

    /**
    * @brief Returns the index of the array element with the smallest absolute value (l1 norm)
    * @param[in] x Input data
    * @return The array index corresponding to the smallest element in the array (0-indexing)
    */
    template<class T> EXPORTCPUCOREMATH size_t amin( hoNDArray< complext<T> > *x );

    /**
    * @brief Returns the index of the array element with the largest absolute value (l1-norm)
    * @param[in] x Input data
    * @return The array index corresponding to the largest element in the array (0-indexing)
    */
    template<class T> EXPORTCPUCOREMATH size_t amax( hoNDArray<T> *x );

    /**
    * @brief Returns the index of the array element with the largest absolute value (l1-norm)
    * @param[in] x Input data
    * @return The array index corresponding to the largest element in the array (0-indexing)
    */
    template<class T> EXPORTCPUCOREMATH size_t amax( hoNDArray< std::complex<T> > *x );

    /**
    * @brief Returns the index of the array element with the largest absolute value (l1-norm)
    * @param[in] x Input data
    * @return The array index corresponding to the largest element in the array (0-indexing)
    */
    template<class T> EXPORTCPUCOREMATH size_t amax( hoNDArray< complext<T> > *x );

    /**
    * @brief Calculates y = a*x+y in which x and y are considered as vectors
    * @param[in] a Scalar value
    * @param[in] x Array
    * @param[in,out] y Array
    */
    template<class T> EXPORTCPUCOREMATH void axpy( T a, hoNDArray<T> *x, hoNDArray<T> *y );
}
