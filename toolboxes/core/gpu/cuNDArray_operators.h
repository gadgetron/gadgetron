/** \file cuNDArray_operators.h
    \brief Common element-wise arithmetic operators on the cuNDArray class.
    
    cuNDArray_operators.h defines element-wise arithmetic array operations on the cuNDArray class.
    We define the common operators +=, -=, *= and \= for both array-array and array-constant operations.
    We have deliberately omitted to define operator+, operator- etc. since this would require returning an cuNDArray,
    in turn invoking an explicit memcpy by the assignment operator.
    Batch mode functionality is provided.
    The implementation is based on Thrust.
    This code is purposely split into a header and underlying implementation (.cu) 
    as this allows specific instantiation of the supported template types. 
    Furthermore thrust code can only be compiled by nvcc.
    The supported types are float, double, Gadgetron::complext<float> and Gadgetron::complext<double>. 
    Scalars can be applied to complex numbers of corresponding precision.
*/

#pragma once

#include "cuNDArray.h"
#include "Gadgetron_enable_types.h"

namespace Gadgetron {

  /**
   * @brief Implementation of element-wise operator+= on two cuNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
   
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<T> & operator+= (cuNDArray<T> &x, const cuNDArray<T> &y);
  
  /**
   * @brief Implementation of element-wise operator+= on a cuNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<T > & operator+= (cuNDArray<T> &x, T y );
    
  /**
   * @brief Implementation of element-wise operator+= on two cuNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
   
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<complext<T > > & operator+= (cuNDArray<complext<T> > &x, const cuNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator+= on a cuNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<complext< T > > & operator+= (cuNDArray<complext<T> > &x, T y );

  /**
   * @brief Implementation of element-wise operator-= on two cuNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
   
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<T > & operator-= (cuNDArray<T> &x, const cuNDArray<T> &y);
  
  /**
   * @brief Implementation of element-wise operator-= on a cuNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<T > & operator-= (cuNDArray<T> &x, T y );
    
  /**
   * @brief Implementation of element-wise operator-= on two cuNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
   
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<complext< T > > & operator-= (cuNDArray<complext<T> > &x, const cuNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator-= on a cuNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<complext< T > > & operator-= (cuNDArray<complext<T> > &x, T y );

  /**
   * @brief Implementation of element-wise operator*= on two cuNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
   
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray< T >  & operator*= (cuNDArray<T> &x, const cuNDArray<T> &y);
  
  /**
   * @brief Implementation of element-wise operator*= on a cuNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T, class = std::enable_if_t<enable_operator<T>()>> cuNDArray<T>  & operator*= (cuNDArray<T> &x, T y );
    
  /**
   * @brief Implementation of element-wise operator*= on two cuNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
   
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<complext< T > > & operator*= (cuNDArray<complext<T> > &x, const cuNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator*= on a cuNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<complext< T > > & operator*= (cuNDArray<complext<T> > &x, T y );

  /**
   * @brief Implementation of element-wise operator/= on two cuNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
   
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray< T >  & operator/= (cuNDArray<T> &x, const cuNDArray<T> &y);
  
  /**
   * @brief Implementation of element-wise operator/= on a cuNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray< T > & operator/= (cuNDArray<T> &x, T y );
    
  /**
   * @brief Implementation of element-wise operator/= on two cuNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
   
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<complext< T > >  & operator/= (cuNDArray<complext<T> > &x, const cuNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator/= on a cuNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T, class = std::enable_if_t<enable_operator_v<T>>> cuNDArray<complext< T > >  & operator/= (cuNDArray<complext<T> > &x, T y );
 /**
   * @brief Implementation of element-wise operator AND on two cuNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.

   * Let y be an n-dimensional array.
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
   cuNDArray<bool >  & operator&= (cuNDArray<bool > &x, cuNDArray<bool> &y);
/**
   * @brief Implementation of element-wise operator OR on two cuNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.

   * Let y be an n-dimensional array.
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  cuNDArray<bool >  & operator|= (cuNDArray<bool > &x, cuNDArray<bool> &y);
}
