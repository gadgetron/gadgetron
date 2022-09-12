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

    /**
     * Creates a 1D histogram of the data. Data outside min_value and max_value will be put in the
     * first and last bin respectively
     * @param data
     * @param bins Number of bins to use
     * @param min_val Minimum value to bin within
     * @param max_val Maximum value to bin within
     * @return
     */
    template <class REAL>
    std::vector<size_t> histogram(const hoNDArray<REAL>& data, size_t bins, REAL min_val,
                                  REAL max_val);
    template <class REAL> std::vector<size_t> histogram(const hoNDArray<REAL>& data, size_t bins);
    /**
     * Calculates a fast, but approximate, percentile of the data.
     * @tparam REAL
     * @param data
     * @param fraction Fraction of the data. Should be in the range [0,1]
     * @param bins Number of bins used for the internal binning
     * @return
     */
    template <class REAL>
    REAL percentile_approx(const hoNDArray<REAL>& data, REAL fraction, size_t bins = 100);

    /**
     * Calculate the percentile of the data, using linear interpolation between points.
     * @tparam REAL
     * @param data
     * @param fraction Fraction of the data. Should be in the range [0,1]
     * @return
     */
    template <class REAL> REAL percentile(const hoNDArray<REAL>& data, REAL fraction);


    /**
     * Calculates the Kullback-Leibler divergence, which measure the difference of distributions between two datasets.
     * @tparam REAL
     * @param dataset1
     * @param dataset2
     * @return
     */
    float jensen_shannon_divergence(const hoNDArray<float>& dataset1, const hoNDArray<float>& dataset2, size_t bins = 100);

    }


#include "hoNDArray_reductions.hxx"