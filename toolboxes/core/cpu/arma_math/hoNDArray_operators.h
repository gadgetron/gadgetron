/** \file hoNDArray_operators.h
    \brief Common element-wise arithmetic operators on the hoNDArray class.
    
    hoNDArray_operators.h defines element-wise arithmetic array operations on the hoNDArray class.
    We define the common operators +=, -=, *= and \= for both array-array and array-constant operations.
    We have deliberately omitted to define operator+, operator- etc. since this would require returning an hoNDArray,
    in turn invoking an explicit memcpy by the assignment operator.
    Batch mode functionality is provided.
    The implementation is based on Armadillo.
    This code is purposely split into a header and underlying implementation (.cpp) 
    as this allows specific instantiation of the supported template types.     
    The supported types are float, double, std::complex<float>, std::complex<double>, 
    Gadgetron::complext<float> and Gadgetron::complext<double>. 
    Scalars can be applied to complex numbers of corresponding precision.
*/

#pragma once

#include "hoNDArray.h"

#include "cpucore_math_export.h"

namespace Gadgetron{

  /**
   * @brief Implementation of element-wise operator+= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray<T>& operator+= (hoNDArray<T> &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator+= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray<T>& operator+= (hoNDArray<T> &x, const T &y);

  /**
   * @brief Implementation of element-wise operator+= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< std::complex<T> >& operator+= (hoNDArray< std::complex<T> > &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator+= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< std::complex<T> >& operator+= (hoNDArray< std::complex<T> >&x, const T &y);

  /**
   * @brief Implementation of element-wise operator+= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< complext<T> >& operator+= (hoNDArray< complext<T> > &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator+= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< complext<T> >& operator+= (hoNDArray< complext<T> >&x, const T &y);

  /**
   * @brief Implementation of element-wise operator-= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray<T>& operator-= (hoNDArray<T> &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator-= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray<T>& operator-= (hoNDArray<T> &x, const T &y);

  /**
   * @brief Implementation of element-wise operator-= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< std::complex<T> >& operator-= (hoNDArray< std::complex<T > > &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator-= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< std::complex<T> >& operator-= (hoNDArray< std::complex<T> >&x, const T &y);

  /**
   * @brief Implementation of element-wise operator-= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< complext<T> >& operator-= (hoNDArray< complext<T > > &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator-= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< complext<T> >& operator-= (hoNDArray< complext<T> >&x, const T &y);

  /**
   * @brief Implementation of element-wise operator*= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray<T>& operator*= (hoNDArray<T> &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator*= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray<T>& operator*= (hoNDArray<T> &x, const T &y);

  /**
   * @brief Implementation of element-wise operator*= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< std::complex<T> >& operator*= (hoNDArray< std::complex<T> > &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator*= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< std::complex<T> >& operator*= (hoNDArray< std::complex<T> > &x, const T &y);

  /**
   * @brief Implementation of element-wise operator*= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< complext<T> >& operator*= (hoNDArray< complext<T> > &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator*= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< complext<T> >& operator*= (hoNDArray< complext<T> > &x, const T &y);

  /**
   * @brief Implementation of element-wise operator/= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray<T>& operator/= (hoNDArray<T> &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator/= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray<T>& operator/= (hoNDArray<T> &x, const T &y);

  /**
   * @brief Implementation of element-wise operator/= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< std::complex<T> >& operator/= (hoNDArray< std::complex<T> > &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator/= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< std::complex<T> >& operator/= (hoNDArray< std::complex<T> > &x, const T &y);

  /**
   * @brief Implementation of element-wise operator/= on two hoNDArrays.
   * @param[in,out] x Input and output array.
   * @param[in] y Input array.
 
   * Let y be an n-dimensional array. 
   * Then the sizes of the first n array dimensions must match between x and y.
   * If x contains further dimensions the operator is batched across those dimensions.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< complext<T> >& operator/= (hoNDArray< complext<T> > &x, const hoNDArray<T> &y);

  /**
   * @brief Implementation of element-wise operator/= on a hoNDArray with a scalar value.
   * @param[in,out] x Input and output array.
   * @param[in] y Input scalar.
   */
  template<class T> EXPORTCPUCOREMATH hoNDArray< complext<T> >& operator/= (hoNDArray< complext<T> > &x, const T &y);
}
