#pragma once

#include "hoNDArray.h"
#include "hoArmadillo.h"
#include "cpucore_math_export.h"

/** \file hoNDArray_ewmath.h
    \brief Element-wise math operations on the hoNDArray class.
    
    hoNDArray_ewmath.h defines element-wise array operations on the hoNDArray class.
    The implementation is based on Armadillo (whenever suitable functions are available).
    Every function comes on two flavours:
    1) A function that returns a smart pointer to a new array holding the result of the element-wise operation, and
    2) A function that perform in-place element-wise computation replacing the input array.
    This code is purposely split into a header and underlying implementation (.cpp) 
    as this allows specific instantiation of the supported template types.     
    The supported types are float, double, std::complex<float>, std::complex<double>, 
    Gadgetron::complext<float> and Gadgetron::complext<double> -- with some deliberate omissions.
*/

namespace Gadgetron{

  /**
   * @brief Calculates the element-wise absolute values (l2 norm) of the array entries
   * @param[in] x Input array
   * @return A new array containing the element-wise absolute values of the input
   */
  template<class T> EXPORTCPUCORE boost::shared_ptr< hoNDArray<typename realType<T>::Type> > abs( hoNDArray<T> *x );

  /**
   * @brief Calculates the element-wise absolute values (l2 norm) of the array entries
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCORE void abs_inplace( hoNDArray<T> *x );
    
  /**
   * @brief Calculates the element-wise sqrt of the array entries
   * @param[in] x Input array
   * @return A new array containing the element-wise sqrt of the input
   */
  template<class T> EXPORTCPUCORE boost::shared_ptr< hoNDArray<T> > sqrt( hoNDArray<T> *x );

  /**
   * @brief Calculates the element-wise sqrt of the array entries
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCORE void sqrt_inplace( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise square of the array entries
   * @param[in] x Input array
   * @return A new array containing the element-wise square of the input
   */
  template<class T> EXPORTCPUCORE boost::shared_ptr< hoNDArray<T> > square( hoNDArray<T> *x );
    
  /**
   * @brief Calculates the element-wise square of the array entries
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCORE void square_inplace( hoNDArray<T> *x );
    
  /**
   * @brief Calculates the element-wise reciprocal of the array entries
   * @param[in] x Input array
   * @return A new array containing the element-wise reciprocal of the input
   */
  template<class T> EXPORTCPUCORE boost::shared_ptr< hoNDArray<T> > reciprocal( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise reciprocal of the array entries
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCORE boost::shared_ptr< hoNDArray<T> > reciprocal_inplace( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise reciprocal sqrt of the array entries
   * @param[in] x Input array
   * @return A new array containing the element-wise reciprocal sqrt of the input
   */
  template<class T> EXPORTCPUCORE boost::shared_ptr< hoNDArray<T> > reciprocal_sqrt( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise reciprocal sqrt of the array entries
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCORE void reciprocal_sqrt_inplace( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the elementwise signum function on the array
   * @param[in] x Input array
   * @return A new array containing the element-wise sgn of the input
   */
  template<class T> EXPORTCPUCORE boost::shared_ptr< hoNDArray<T> > sgn( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the elementwise signum function on the array
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCORE void sgn_inplace( hoNDArray<T> *x );

  /**
   * @brief Clamps all values in the array to the minimum and maximum values specified
   * @param[in] x Input array
   * @param[in] min minimum value
   * @param[in] max maximum value
   * @return A new array containing the element-wise clamped values of the input
   */
  template<class T> EXPORTCPUCORE boost::shared_ptr< hoNDArray<T> > clamp( hoNDArray<T> *x, T min, T max );

  /**
   * @brief Clamps all values in the array to the minimum and maximum values specified
   * @param[in,out] x Input and output array
   * @param[in] min minimum value
   * @param[in] max maximum value
   */
  template<class T> EXPORTCPUCORE void clamp_inplace( hoNDArray<T> *x, T min, T max );

  /**
   * @brief Clamps all values in the array to a minimum value allowed.
   * @param[in] x Input array
   * @param[in] min Minimum value
   * @return A new array containing the element-wise clamped values of the input
   */
  template<class T> EXPORTCPUCORE boost::shared_ptr< hoNDArray<T> > clamp_min( hoNDArray<T> *x, T min );

  /**
   * @brief Clamps all values in the array to a minimum value allowed.
   * @param[in,out] x Input and output array
   * @param[in] min Minimum value
   */
  template<class T> EXPORTCPUCORE void clamp_min_inplace( hoNDArray<T> *x, T min );

  /**
   * @brief Clamps all values in the array to a maximum value allowed.
   * @param[in] x Input array
   * @param[in] max Maximum value
   * @return A new array containing the element-wise clamped values of the input
   */
  template<class T> EXPORTCPUCORE boost::shared_ptr< hoNDArray<T> > clamp_max( hoNDArray<T> *x, T max );

  /**
   * @brief Clamps all values in the array to a maximum value allowed.
   * @param[in,out] x Input and output array
   * @param[in] max Maximum value
   */
  template<class T> EXPORTCPUCORE void clamp_max_inplace( hoNDArray<T> *x, T max );
}
