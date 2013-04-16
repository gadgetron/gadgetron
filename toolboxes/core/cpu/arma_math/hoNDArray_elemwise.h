/** \file hoNDArray_elemwise.h
    \brief Element-wise math operations on the hoNDArray class.
    
    hoNDArray_elementwise.h defines element-wise array operations on the hoNDArray class.
    Many of the provided functions come in two flavours:
    1) A function that returns a smart pointer to a new array holding the result of the element-wise operation, and
    2) A function that perform in-place element-wise computation replacing the input array.
    When both versions are available the in-place version is suffixed _inplace.
    Some functions (clear, fill, clamp, clamp_min, clamp_max, normalize, shrink1, shrinkd) are only provided as in-place operations,
    and they do not carry the _inplace suffix in order to keep user code compact.
    A few functions return a different type as its input array 
    (abs on complex data, real, imag, real_to_std_complex, real_to_complext) and consequently is not offered as an in place operation.
    The functions provided in hoNDArray_elemwise are deliberatly placed outside the NDArray derived classes
    - to allow the NDArray classes to be lightweight header only data containers for both the cpu and gpu instances
    - to allow for external library optimized implementations of the element-wise functions without adding such dependencies to the core data container
    The present cpu implementation is based on Armadillo (whenever suitable functions are available).
    The implementation is purposely split into a header and underlying implementation (.cpp) 
    as this allows specific instantiation of the supported template types.     
    The supported types are float, double, std::complex<float>, std::complex<double>, 
    Gadgetron::complext<float> and Gadgetron::complext<double> -- with some deliberate omissions.
*/

#pragma once

#include "hoNDArray.h"
#include "cpucore_math_export.h"

namespace Gadgetron{

  /**
   * @brief Calculates the element-wise absolute values (l2 norm) of the array entries
   * @param[in] x Input array.
   * @return A new array containing the element-wise absolute values of the input.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<typename realType<T>::Type> > abs( hoNDArray<T> *x );

  /**
   * @brief Calculates the element-wise absolute values (l2 norm) of the array entries (in place).
   * @param[in,out] x Input and output array.
   */
  template<class T> EXPORTCPUCOREMATH void abs_inplace( hoNDArray<T> *x );

  /**
   * @brief Calculates the element-wise squared absolute values of the array entries
   * @param[in] x Input array.
   * @return A new array containing the element-wise absolute values of the input.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<typename realType<T>::Type> > abs_square( hoNDArray<T> *x );
    
  /**
   * @brief Calculates the element-wise sqrt of the array entries.
   * @param[in] x Input array.
   * @return A new array containing the element-wise sqrt of the input.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > sqrt( hoNDArray<T> *x );

  /**
   * @brief Calculates the element-wise sqrt of the array entries (in place).
   * @param[in,out] x Input and output array.
   */
  template<class T> EXPORTCPUCOREMATH void sqrt_inplace( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise square of the array entries.
   * @param[in] x Input array.
   * @return A new array containing the element-wise square of the input.
   *
   * For real numbers this functions is equivalent to square. 
   * For complex arrays abs_square() and square() differ however.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > square( hoNDArray<T> *x );
    
  /**
   * @brief Calculates the element-wise square of the array entries (in place).
   * @param[in,out] x Input and output array.
   */
  template<class T> EXPORTCPUCOREMATH void square_inplace( hoNDArray<T> *x );
    
  /**
   * @brief Calculates the element-wise reciprocal of the array entries.
   * @param[in] x Input array.
   * @return A new array containing the element-wise reciprocal of the input.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > reciprocal( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise reciprocal of the array entries (in place).
   * @param[in,out] x Input and output array.
   */
  template<class T> EXPORTCPUCOREMATH void reciprocal_inplace( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise reciprocal sqrt of the array entries.
   * @param[in] x Input array.
   * @return A new array containing the element-wise reciprocal sqrt of the input.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > reciprocal_sqrt( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the element-wise reciprocal sqrt of the array entries (in place).
   * @param[in,out] x Input and output array.
   */
  template<class T> EXPORTCPUCOREMATH void reciprocal_sqrt_inplace( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the elementwise signum function on the array.
   * @param[in] x Input array.
   * @return A new array containing the element-wise sgn of the input.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > sgn( hoNDArray<T> *x );
  
  /**
   * @brief Calculates the elementwise signum function on the array (in place).
   * @param[in,out] x Input and output array.
   */
  template<class T> EXPORTCPUCOREMATH void sgn_inplace( hoNDArray<T> *x );

  /**
   * @brief Extract the real component from a complex array.
   * @param[in] x Input array.
   * @return A new array of the real component of the complex array.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<typename realType<T>::Type> > real( hoNDArray<T> *x );

  /**
   * @brief Extract the imaginary component from a complex array.
   * @param[in] x Input array.
   * @return A new array of the imaginary component of the complex array.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<typename realType<T>::Type> > imag( hoNDArray<T> *x );

  /**
   * @brief Create a new array of the complex conjugate of the input array. For real arrays a copy of the input array is return.
   * @param[in] x Input array.
   * @return A new array of the complex conjugate of the input array.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > conj( hoNDArray<T> *x );

  /**
   * @brief Construct a complex array from a real array.
   * @param[in] x Input array.
   * @return A new complex array containing the input array in the real component and zeros in the imaginary component.
   */
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > real_to_complex( hoNDArray<typename realType<T>::Type> *x );
  
  //
  // From hereon the functions are all in-place although without the _inplace suffix...
  //

  /**
   * @brief Clears the array to all zeros ( in place). Faster than fill.
   * @param[in,out] x Input and output array.
   */
  template<class T> EXPORTCPUCOREMATH void clear( hoNDArray<T> *x );

  /**
   * @brief Fills the array with a user provided constant value (in place).
   * @param[in,out] x Input and output array.
   * @param[in] val Fill value.
   */
  template<class T> EXPORTCPUCOREMATH void fill( hoNDArray<T> *x, T val );

  /**
   * @brief Clamps all values in the array to the minimum and maximum values specified (in place).
   * @param[in,out] x Input and output array.
   * @param[in] min minimum value.
   * @param[in] max maximum value.
   */
  template<class T> EXPORTCPUCOREMATH void clamp( hoNDArray<T> *x, typename realType<T>::Type min, typename realType<T>::Type max );
  
  /**
   * @brief Clamps all values in the array to a minimum value allowed (in place).
   * @param[in,out] x Input and output array.
   * @param[in] min Minimum value.
   */
  template<class T> EXPORTCPUCOREMATH void clamp_min( hoNDArray<T> *x, typename realType<T>::Type min );

  /**
   * @brief Clamps all values in the array to a maximum value allowed (in place).
   * @param[in,out] x Input and output array.
   * @param[in] max Maximum value.
   */
  template<class T> EXPORTCPUCOREMATH void clamp_max( hoNDArray<T> *x, typename realType<T>::Type max );

  /**
   * @brief In place normalization (scaling) to a new maximum absolute array value val.
   * @param[in,out] x Input and output array.
   * @param[in] val New maximum absolute array value (according to the l2-norm)
   */  
  template<class T> EXPORTCPUCOREMATH void normalize( hoNDArray<T> *x, typename realType<T>::Type val = typename realType<T>::Type(1) );

  /**
   * @brief Shrinkage (soft thresholding), i.e. shrink(x,gamma) = x/abs(x)*max(abs(x)-gamma,0).
   * @param[out] out Output array. Can be 0x0 in which case an in place transform is performed.
   * @param[in,out] x Input array (and output array if out == 0x0).
   * @param[in] gamma Shrinkage control parameter
   */  
  template<class T> EXPORTCPUCOREMATH void shrink1( hoNDArray<T> *x, typename realType<T>::Type gamma, hoNDArray<T> *out = 0x0 );

  /**
   * @brief Shrinkage (soft thresholding, multi-dimensional), i.e. shrink(x,gamma,s) = x/s*max(s-gamma,0).
   * @param[out] out Output array. Can be 0x0 in which case an in place transform is performed.
   * @param[in,out] x Input array (and output array if out == 0x0).
   * @param[in] s Input array, normalization.
   * @param[in] gamma Shrinkage control parameter
   */  
  template<class T> EXPORTCPUCOREMATH void shrinkd ( hoNDArray<T> *x, hoNDArray<typename realType<T>::Type> *s, typename realType<T>::Type gamma, hoNDArray<T> *out = 0x0 );
}
