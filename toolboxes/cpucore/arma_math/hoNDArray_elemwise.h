/** \file hoNDArray_elemwise.h
    \brief Element-wise math operations on the hoNDArray class.
    
    hoNDArray_elementwise.h defines element-wise array operations on the hoNDArray class.
    Most functions come in two flavours:
    1) A function that returns a smart pointer to a new array holding the result of the element-wise operation, and
    2) A function that perform in-place element-wise computation replacing the input array.
    When both versions are available the in-place version is suffixed _inplace.
    A few functions (clear, fill, clamp, clamp_min, clamp_max) are only provided as in-place operations,
    and they do not carry the _inplace suffix in order to keep user code compact.
    The functions provided in hoNDArray_elemwise are deliberatly placed outside the NDArray derived class
    - to allow the NDArray classes to be lightweight header only data containers for both the cpu and gpu instances
    - to allow for external library optimized implementations without adding such dependencies to the core data container
    The present cpu implementation is based on Armadillo (whenever suitable functions are available).
    The implementation is purposely split into a header and underlying implementation (.cpp) 
    as this allows specific instantiation of the supported template types.     
    The supported types are float, double, std::complex<float>, std::complex<double>, 
    Gadgetron::complext<float> and Gadgetron::complext<double> -- with some deliberate omissions.
*/

#pragma once

#include "hoNDArray.h"
#include "hoArmadillo.h"
#include "cpucore_math_export.h"

namespace Gadgetron{

  /**
   * @brief Calculates the element-wise absolute values (l2 norm) of the array entries
   * @param[in] x Input array
   * @return A new array containing the element-wise absolute values of the input
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<typename realType<T>::Type> > abs( hoNDArray<T> *x );

  /**
   * @brief Calculates the element-wise absolute values (l2 norm) of the array entries
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCOREMATH void abs_inplace( hoNDArray<T> *x );
    
  /**
   * @brief Calculates the element-wise sqrt of the array entries
   * @param[in] x Input array
   * @return A new array containing the element-wise sqrt of the input
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > sqrt( hoNDArray<T> *x );

  /**
   * @brief Calculates the element-wise sqrt of the array entries
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCOREMATH void sqrt_inplace( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise square of the array entries
   * @param[in] x Input array
   * @return A new array containing the element-wise square of the input
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > square( hoNDArray<T> *x );
    
  /**
   * @brief Calculates the element-wise square of the array entries
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCOREMATH void square_inplace( hoNDArray<T> *x );
    
  /**
   * @brief Calculates the element-wise reciprocal of the array entries
   * @param[in] x Input array
   * @return A new array containing the element-wise reciprocal of the input
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > reciprocal( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise reciprocal of the array entries
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCOREMATH void reciprocal_inplace( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise reciprocal sqrt of the array entries
   * @param[in] x Input array
   * @return A new array containing the element-wise reciprocal sqrt of the input
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > reciprocal_sqrt( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise reciprocal sqrt of the array entries
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCOREMATH void reciprocal_sqrt_inplace( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the elementwise signum function on the array
   * @param[in] x Input array
   * @return A new array containing the element-wise sgn of the input
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > sgn( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the elementwise signum function on the array
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCOREMATH void sgn_inplace( hoNDArray<T> *x );

  //
  // From hereon the functions are all in-place although without the _inplace suffix...
  //

  /**
   * @brief Clears the array to all zeros. Faster than fill.
   * @param[in,out] x Input and output array
   */
  template<class T> EXPORTCPUCOREMATH void clear( hoNDArray<T> *x );

  /**
   * @brief Fills the array with a user provided constant value
   * @param[in,out] x Input and output array
   * @param[in] val Fill value
   */
  template<class T> EXPORTCPUCOREMATH void fill( hoNDArray<T> *x, T val );

  /**
   * @brief Clamps all values in the array to the minimum and maximum values specified
   * @param[in,out] x Input and output array
   * @param[in] min minimum value
   * @param[in] max maximum value
   */
  template<class T> EXPORTCPUCOREMATH void clamp( hoNDArray<T> *x, T min, T max );

  /**
   * @brief Clamps all values in the array to a minimum value allowed.
   * @param[in,out] x Input and output array
   * @param[in] min Minimum value
   */
  template<class T> EXPORTCPUCOREMATH void clamp_min( hoNDArray<T> *x, T min );

  /**
   * @brief Clamps all values in the array to a maximum value allowed.
   * @param[in,out] x Input and output array
   * @param[in] max Maximum value
   */
  template<class T> EXPORTCPUCOREMATH void clamp_max( hoNDArray<T> *x, T max );
}
