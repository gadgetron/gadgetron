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

    /***
    * Calculates the variance of the array
    */
    template<class T> EXPORTCPUCOREMATH T var(hoNDArray<T>* data);

    /***
     * Calulates the median of the array
     */
    template<class T> EXPORTCPUCOREMATH T median(hoNDArray<T>* data);

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
    template<class T> EXPORTCPUCOREMATH void asum(const hoNDArray<T>& x, typename realType<T>::Type& r);
    template<class T> EXPORTCPUCOREMATH typename realType<T>::Type asum(const hoNDArray<T>& x);

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
    * @brief ind = min(abs(x(:))
    find the minimal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T> EXPORTCPUCOREMATH 
    void minAbsolute(const hoNDArray<T>& x, T& r, size_t& ind);

    /**
    * @brief ind = max(abs(x(:))
    find the miximal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T> EXPORTCPUCOREMATH 
    void maxAbsolute(const hoNDArray<T>& x, T& r, size_t& ind);

    /**
    * @brief r = norm(x(:), 2)
    compute L2 norm of x
    */
    template <typename T> EXPORTCPUCOREMATH 
    void norm2(const hoNDArray<T>& x, typename realType<T>::Type& r);

    template <typename T> EXPORTCPUCOREMATH 
    typename realType<T>::Type norm2(const hoNDArray<T>& x);

    /**
    * @brief r = norm(x(:), 1)
    compute L1 norm of x = sum( abs(x(:) )
    */
    template <typename T> EXPORTCPUCOREMATH 
    void norm1(const hoNDArray<T>& x, typename realType<T>::Type& r);

    template <typename T> EXPORTCPUCOREMATH 
    typename realType<T>::Type norm1(const hoNDArray<T>& x);

    /**
    * @brief dot product of x and conj(y)
    r = x dot conj(y)
    */
    template <typename T> EXPORTCPUCOREMATH 
    void dotc(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r);

    template <typename T> EXPORTCPUCOREMATH 
    T dotc(const hoNDArray<T>& x, const hoNDArray<T>& y);

    /**
    * @brief dot product of x and y
    r = x dot y
    */
    template <typename T> EXPORTCPUCOREMATH 
    void dotu(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r);

    template <typename T> EXPORTCPUCOREMATH 
    T dotu(const hoNDArray<T>& x, const hoNDArray<T>& y);

    /**
    * @brief sort the ND array
    */
    template <typename T> EXPORTCPUCOREMATH void sort(const hoNDArray<T>& x, hoNDArray<T>& r, bool isascending);

    /**
    * @brief sort the ND array with indexes, similar to [r, ind] = sor(x(:))
    */
    template <typename T> EXPORTCPUCOREMATH void sort(const hoNDArray<T>& x, hoNDArray<T>& r, std::vector<size_t>& ind, bool isascending);

    /**
    * @brief finds the index of the element with the maximal absolute value.
    */
    template<class T> EXPORTCPUCOREMATH size_t amax(const hoNDArray<T>& x);

    /**
    * @brief get the min and max value from an array (only for float and double type)
    */
    template <class T> EXPORTCPUCOREMATH 
    void minValue(const hoNDArray<T>& a, T& v);

    template <class T> EXPORTCPUCOREMATH 
    void maxValue(const hoNDArray<T>& a, T& v);
}
