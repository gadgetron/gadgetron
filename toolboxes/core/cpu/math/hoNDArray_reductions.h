#pragma once

#include "hoNDArray.h"
#include "cpp_blas.h"

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
    template<class REAL>  REAL max(const hoNDArray<REAL>* data);
    template<class REAL>  REAL max(const hoNDArray<REAL>& data);

    /***
    * Finds the minimum element of the array
    */
    template<class REAL>  REAL min(const hoNDArray<REAL>* data);
    template<class REAL>  REAL min(const hoNDArray<REAL>& data);

    /***
    * Finds the mean of the array
    */
    template<class T>  T mean(const hoNDArray<T>* data);
    template<class T>  T mean(const hoNDArray<T>& data);

    /***
    * Calculates the sum of the array
    */
    template<class T>  T sum(const hoNDArray<T>* data);
    template<class T>  T sum(const hoNDArray<T>& data);

    /***
    * Calculates the std of the array
    */
    template<class T>  T stddev(const hoNDArray<T>* data);
    template<class T>  T stddev(const hoNDArray<T>& data);

    /***
    * Calculates the variance of the array
    */
    template<class T>  T var(const hoNDArray<T>* data);
    template<class T>  T var(const hoNDArray<T>& data);

    /***
     * Calulates the median of the array
     */
    template<class T>  T median(const hoNDArray<T>* data);
    template<class T>  T median(const hoNDArray<T>& data);

    /**
    * @brief Calculates the dot product of two arrays (as vectors).
    * @param[in] x Array 1. For complex arrays the complex conjugate of x is used.
    * @param[in] y Array 2.
    * @param[in] cc Specifies whether to use the complex conjugate of x (when applicable).
    * @return The dot product of x and y
    */
    template<class T> T dot(const hoNDArray<T> *x, const hoNDArray<T> *y, bool cc = true );
    template<class T> T dot(const hoNDArray<T>& x, const hoNDArray<T>& y, bool cc = true );


    /**
    * @brief Calculates the sum of the l1-norms of the array entries
    * @param[in] arr Input array
    * @return The l1-norm of the array
    */
    template<class T> typename realType<T>::Type asum( const hoNDArray<T> *x );
    template<class T> typename realType<T>::Type asum(const hoNDArray<T>& x);

    /**
    * @brief Calculates the l2-norm of the array (as a vector)
    * @param[in] arr Input array
    * @return The l2-norm of the array
    */
    template<class T> typename realType<T>::Type nrm2(const hoNDArray<T> *x );
    template<class T> typename realType<T>::Type nrm2(const hoNDArray<T>& x );


    /**
    * @brief Returns the index of the array element with the smallest absolute value (l1 norm)
    * @param[in] x Input data
    * @return The array index corresponding to the smallest element in the array (0-indexing)
    */
    template<class T>  size_t amin( const hoNDArray<T> *x );



    /**
    * @brief Returns the index of the array element with the largest absolute value (l1-norm)
    * @param[in] x Input data
    * @return The array index corresponding to the largest element in the array (0-indexing)
    */
    template<class T> size_t amax(const  hoNDArray<T> *x );

   /**
    * @brief finds the index of the element with the maximal absolute value.
    */
    template<class T> size_t amax(const hoNDArray<T>& x);

    /**
    * @brief ind = min(abs(x(:))
    find the minimal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T>  
    void minAbsolute(const hoNDArray<T>& x, T& r, size_t& ind);

    /**
    * @brief ind = max(abs(x(:))
    find the miximal absolute value of x and its position index ind
    r = x[ind], not abs(x[ind])
    */
    template <typename T>  
    void maxAbsolute(const hoNDArray<T>& x, T& r, size_t& ind);





    /**
    * @brief sort the ND array
    */
    template <typename T>  void sort(const hoNDArray<T>& x, hoNDArray<T>& r, bool isascending);

    /**
    * @brief sort the ND array with indexes, similar to [r, ind] = sor(x(:))
    */
    template <typename T>  void sort(const hoNDArray<T>& x, hoNDArray<T>& r, std::vector<size_t>& ind, bool isascending);



    /**
    * @brief get the min and max value from an array (only for float and double type)
    */
    template <class T>  
    void minValue(const hoNDArray<T>& a, T& v);

    template <class T>  
    void maxValue(const hoNDArray<T>& a, T& v);


    template<class REAL>
    REAL percentile_approx(const hoNDArray<REAL>& data, REAL fraction,size_t bins = 100);

     template<class REAL>
    REAL percentile(const hoNDArray<REAL>& data, REAL fraction);

}


#include "hoNDArray_reductions.hxx"