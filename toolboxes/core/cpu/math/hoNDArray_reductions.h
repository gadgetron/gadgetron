#pragma once

#include "hoNDArray.h"
#include "cpucore_math_export.h"

#ifdef max
    #undef max
#endif // max

#ifdef min
    #undef min
#endif // min

namespace Gadgetron{

    /***
    * Finds the maximum element of the array
    */
    template<class REAL> EXPORTCPUCOREMATH REAL max(hoNDArray<REAL>* data);

    /***
    * Finds the minimum element of the array
    */
    template<class REAL> EXPORTCPUCOREMATH REAL min(hoNDArray<REAL>* data);

    /***
    * Finds the mean of the array
    */
    template<class T> EXPORTCPUCOREMATH T mean(hoNDArray<T>* data);

    /***
    * Calculates the sum of the array
    */
    template<class T> EXPORTCPUCOREMATH T sum(hoNDArray<T>* data);

    /***
    * Calculates the std of the array
    */
    template<class T> EXPORTCPUCOREMATH T stddev(hoNDArray<T>* data);

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
}
