#pragma once
#include "complext.h"
#include "hoNDArray.h"
namespace Gadgetron{
//TODO: Reimplement using BLAS

/**
 * @brief Calculates the dot product of two arrays
 * @param[in] x Array 1
 * @param[in] y Array 2
 * @return The dot product of x and y
 */
template<class T> T dot(hoNDArray<T> *x,hoNDArray<T> *y);

/**
 * @brief Calculates the l2-norm of the array
 * @param[in] arr Input array
 * @return
 */
template<class T> typename realType<T>::type nrm2( hoNDArray<T>* arr);

/**
 * @brief Calculates Y = a*X+Y
 * @param[in] a Scalar
 * @param[in] x Array
 * @param[in,out] y Array
 */
template<class T> void axpy(T a, hoNDArray<T>* x, hoNDArray<T>* y);

template<class T> void axpy(T a, hoNDArray<complext<T> >* x, hoNDArray<complext<T> >* y){
	axpy(complext<T>(a),x,y);
}

/**
 * @brief Gets the index of the index of the element with minimum absolute
 * @param[in] x Input data
 * @return index of absolute minimum values
 */
template<class T> int amin(hoNDArray<T>* x);
/**
 * @brief Gets the index of the index of the element with minimum absolute
 * @param x Input data
 * @return index of absolute minimum values
 */
template<class T> int amin(hoNDArray<T>* x);

/**
 * @brief Gets the index of the index of the element with maximum absolute
 * @param[in] x Input data
 * @return index of absolute maximum values
 */
template<class T> int amax(hoNDArray<T>* x);


template<class T> typename realType<T>::type asum(hoNDArray<T>* x);

}
