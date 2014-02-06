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

#include "GadgetronCommon.h"
#include <complex>

#ifdef USE_MKL
#include "mkl.h"
#endif // USE_MKL

#ifdef GT_Complex8
#undef GT_Complex8
#endif // GT_Complex8
typedef std::complex<float> GT_Complex8;

#ifdef GT_Complex16
#undef GT_Complex16
#endif // GT_Complex16
typedef std::complex<double> GT_Complex16;

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
  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > 
  real_to_complex( hoNDArray<typename realType<T>::Type> *x );

  template<class T> EXPORTCPUCOREMATH boost::shared_ptr< hoNDArray<T> > 
  real_imag_to_complex( hoNDArray<typename realType<T>::Type> *real, hoNDArray<typename realType<T>::Type>* imag);

  template<class T> EXPORTCPUCOREMATH bool 
  real_imag_to_complex(const hoNDArray<typename realType<T>::Type>& real, const hoNDArray<typename realType<T>::Type>& imag, hoNDArray<T>& cplx);

  template<class T> EXPORTCPUCOREMATH bool 
  complex_to_real_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real, hoNDArray<typename realType<T>::Type>& imag);

  template<> EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray<float>& cplx, hoNDArray<float>& real, hoNDArray<float>& imag);
  template<> EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray<double>& cplx, hoNDArray<double>& real, hoNDArray<double>& imag);

  template<class T> EXPORTCPUCOREMATH bool complex_to_real(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real);
  template<class T> EXPORTCPUCOREMATH bool complex_to_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& imag);

  //
  // From hereon the functions are all in-place although without the _inplace suffix...
  //

  /**
   * @brief Clears the array to all zeros ( in place). Faster than fill.
   * @param[in,out] x Input and output array.
   */
  template<class T> EXPORTCPUCOREMATH void clear( hoNDArray<T> *x );
  template<class T> EXPORTCPUCOREMATH void clear( hoNDArray<T>& x );

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
   * @param[in] min_val value to which everything below the minimum will be set
   * @param[in] max_val value to which everything above the maximum will be set
   */
  template<class T> EXPORTCPUCOREMATH void clamp( hoNDArray<T> *x, typename realType<T>::Type min, typename realType<T>::Type max, T min_val, T max_val );
  
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
   * @brief In place p-shrinkage (soft thresholding), i.e. pshrink(x,gamma,p) = x/abs(x)*max(abs(x)-gamma*abs(x)^(p-1),0).
   * @param[out] out Output array. Can be 0x0 in which case an in place transform is performed.
   * @param[in,out] x Input array (and output array if out == 0x0).
   * @param[in] gamma Shrinkage control parameter
   * @param[in] p p value of the shrinkage. Should be less than 1 and more than 0.
   */
	template<class T> EXPORTCPUCOREMATH void pshrink( hoNDArray<T> *x, typename realType<T>::Type gamma,typename realType<T>::Type p, hoNDArray<T> *out = 0x0 );

  /**
   * @brief Shrinkage (soft thresholding, multi-dimensional), i.e. shrink(x,gamma,s) = x/s*max(s-gamma,0).
   * @param[out] out Output array. Can be 0x0 in which case an in place transform is performed.
   * @param[in,out] x Input array (and output array if out == 0x0).
   * @param[in] s Input array, normalization.
   * @param[in] gamma Shrinkage control parameter
   */  
  template<class T> EXPORTCPUCOREMATH void shrinkd ( hoNDArray<T> *x, hoNDArray<typename realType<T>::Type> *s, typename realType<T>::Type gamma, hoNDArray<T> *out = 0x0 );

  /**
   * @brief In place p-shrinkage (soft thresholding, multi-dimensional), i.e. pshrink(x,s,gamma,p) = x/s*max(s-gamma*s^(p-1),0).
   * @param[out] out Output array. Can be 0x0 in which case an in place transform is performed.
   * @param[in,out] x Input array (and output array if out == 0x0).
   * @param[in] gamma Shrinkage control parameter
   * @param[in] p p value of the shrinkage. Should be less than 1 and more than 0.
   */
  template<class T> EXPORTCPUCOREMATH void pshrinkd ( hoNDArray<T> *x, hoNDArray<typename realType<T>::Type> *s, typename realType<T>::Type gamma,typename realType<T>::Type p, hoNDArray<T> *out = 0x0 );

#ifdef USE_MKL

    // besides the arma calls, some functions are implemented with the MKL vector utilities

    EXPORTCPUCOREMATH bool add(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r); // r = x + y
    EXPORTCPUCOREMATH bool subtract(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r); // r = x - y
    EXPORTCPUCOREMATH bool multiply(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r); // r = x * y
    EXPORTCPUCOREMATH bool divide(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r); // r = x / y
    EXPORTCPUCOREMATH bool absolute(const hoNDArray<float>& x, hoNDArray<float>& r); // r = abs(x)
    EXPORTCPUCOREMATH bool argument(const hoNDArray<float>& x, hoNDArray<float>& r); // r = angle(x)
    EXPORTCPUCOREMATH bool sqrt(const hoNDArray<float>& x, hoNDArray<float>& r); // r = sqrt(x)
    EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<float>& x, float& r, size_t& ind); // minimal absolute value and index
    EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<float>& x, float& r, size_t& ind); // maximal absolute value and index
    EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<float>& x); // x = x + Epsilon if x==0, prepare for division
    EXPORTCPUCOREMATH bool norm2(const hoNDArray<float>& x, float& r);
    EXPORTCPUCOREMATH bool norm1(const hoNDArray<float>& x, float& r);
    EXPORTCPUCOREMATH bool conv2(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& z); // x: input data, y: convolution kernel, z: output; each 2D slice is convolved
    EXPORTCPUCOREMATH bool conv3(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& z); // x: input data, y: convolution kernel, z: output; each 3D volume is convolved
    EXPORTCPUCOREMATH bool inv(const hoNDArray<float>& x, hoNDArray<float>& r); // r = 1/x

    EXPORTCPUCOREMATH bool add(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool subtract(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool multiply(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool divide(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool absolute(const hoNDArray<double>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool argument(const hoNDArray<double>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool sqrt(const hoNDArray<double>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<double>& x, double& r, size_t& ind);
    EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<double>& x, double& r, size_t& ind);
    EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<double>& x);
    EXPORTCPUCOREMATH bool norm2(const hoNDArray<double>& x, double& r);
    EXPORTCPUCOREMATH bool norm1(const hoNDArray<double>& x, double& r);
    EXPORTCPUCOREMATH bool conv2(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& z);
    EXPORTCPUCOREMATH bool conv3(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& z);
    EXPORTCPUCOREMATH bool inv(const hoNDArray<double>& x, hoNDArray<double>& r);

    EXPORTCPUCOREMATH bool add(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool subtract(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool multiply(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool divide(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r);
    EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool sqrt(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind);
    EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind);
    EXPORTCPUCOREMATH bool multiplyConj(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r); // r = x * conj(y)
    EXPORTCPUCOREMATH bool argument(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r); // r = angle(x)
    EXPORTCPUCOREMATH bool conjugate(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r); // r = conj(x)
    EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<GT_Complex8>& x);
    EXPORTCPUCOREMATH bool norm2(const hoNDArray<GT_Complex8>& x, float& r);
    EXPORTCPUCOREMATH bool norm1(const hoNDArray<GT_Complex8>& x, float& r);
    EXPORTCPUCOREMATH bool dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, GT_Complex8& r); // x'*y, x and y are N*1 vector
    EXPORTCPUCOREMATH bool conv2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& z);
    EXPORTCPUCOREMATH bool conv3(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& z);
    EXPORTCPUCOREMATH bool inv(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);

    EXPORTCPUCOREMATH bool add(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool subtract(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool multiply(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool divide(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool sqrt(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind);
    EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind);
    EXPORTCPUCOREMATH bool multiplyConj(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool argument(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r);
    EXPORTCPUCOREMATH bool conjugate(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);
    EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<GT_Complex16>& x);
    EXPORTCPUCOREMATH bool norm2(const hoNDArray<GT_Complex16>& x, double& r);
    EXPORTCPUCOREMATH bool norm1(const hoNDArray<GT_Complex16>& x, double& r);
    EXPORTCPUCOREMATH bool dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, GT_Complex16& r);
    EXPORTCPUCOREMATH bool conv2(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& z);
    EXPORTCPUCOREMATH bool conv3(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& z);
    EXPORTCPUCOREMATH bool inv(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template<typename T> EXPORTCPUCOREMATH bool sumOverLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // e.g. for a 4D array, sum over the 4th dimension
    template<typename T> EXPORTCPUCOREMATH bool sumOverSecondLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // e.g. for a 4D array, sum over the 3rd dimension

    template<typename T> EXPORTCPUCOREMATH bool multiplyOverLastDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r); // e.g. x is 3D and y is 4D array, r(:,:,:,n) = y(:,:,:,n) .* x
    template<typename T> EXPORTCPUCOREMATH bool divideOverLastDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r); // e.g. x is 3D and y is 4D array, r(:,:,:,n) = y(:,:,:,n) ./ x

    template<typename T> EXPORTCPUCOREMATH bool sumOver1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // e.g. for a 2D array, sum over the 1st dimension and get an array of [1 E1]
    template<typename T> EXPORTCPUCOREMATH bool sumOver2ndDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // e.g. for a 3D array, sum over the 2nd dimension and get an array of [RO 1 CHA]
    template<typename T> EXPORTCPUCOREMATH bool sumOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // e.g. for a 4D array, sum over the 3rd dimension and get an array of [RO E1 1 N]
    template<typename T> EXPORTCPUCOREMATH bool sumOver4thDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // e.g. for a 5D array [RO E1 CHA N S], sum over the 4th dimension and get an array of [RO E1 CHA 1 S]
    template<typename T> EXPORTCPUCOREMATH bool sumOver5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // e.g. for a 6D array, sum over the 5th dimension and get an array [RO E1 CHA N 1 P]

    template<typename T> EXPORTCPUCOREMATH bool multiplyOver3rdDimension(const hoNDArray<T>& x3D, const hoNDArray<T>& y4D, hoNDArray<T>& r); // e.g. x is 3D and y is 4D array, r(:,:,n,:) = y(:,:,n,:) .* x
    template<typename T> EXPORTCPUCOREMATH bool multiplyOver4thDimension(const hoNDArray<T>& x4D, const hoNDArray<T>& y5D, hoNDArray<T>& r); // e.g. x is 4D and y is 5D array, r(:,:,:,n,:) = y(:,:,:,n,:) .* x
    template<typename T> EXPORTCPUCOREMATH bool multiplyOver5thDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r); // e.g. x is 5D and y is 6D array, r(:,:,:,:, n,:) = y(:,:,:,:,n,:) .* x

    template<typename T> EXPORTCPUCOREMATH bool multiplyOver4thDimensionExcept(const hoNDArray<T>& x4D, const hoNDArray<T>& y5D, size_t n, hoNDArray<T>& r, bool copyY2R=true); // e.g. x is 4D and y is 5D array, r(:,:,:,t,:) = y(:,:,:,t,:) .* x, except for r(:,:,:,n,:) = y(:,:,:,n,:)
    template<typename T> EXPORTCPUCOREMATH bool multiplyOver5thDimensionExcept(const hoNDArray<T>& x, const hoNDArray<T>& y, size_t n, hoNDArray<T>& r, bool copyY2R=true); // e.g. x is 5D and y is 6D array, r(:,:,:,:,t,:) = y(:,:,:,:,t,:) .* x, except for r(:,:,:,:,n,:) = y(:,:,:,:,n,:)

    template<typename T> EXPORTCPUCOREMATH bool multipleAdd(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r); // r = x + y for every part of y
    template<typename T> EXPORTCPUCOREMATH bool multipleMultiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r); // r = x * y for every part of y
    template<typename T> EXPORTCPUCOREMATH bool multipleDivide(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r); // r = x / y for every part of y

    template<typename T> EXPORTCPUCOREMATH bool cropUpTo10DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const std::vector<size_t>& startND, std::vector<size_t>& size);
    template<typename T> EXPORTCPUCOREMATH bool setSubArrayUpTo10DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const std::vector<size_t>& startND, std::vector<size_t>& size);

    template<typename T> EXPORTCPUCOREMATH bool cropOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end);
    template<typename T> EXPORTCPUCOREMATH bool setSubArrayOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end);

    template<typename T> EXPORTCPUCOREMATH bool stdOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& std, bool NMinusOne); // compute the standard deviation along the 3rd dimension, if NMinusOne == true, divided by N-1; otherwise, divided by N

    // template<typename T> EXPORTCPUCOREMATH bool permuteLastTwoDimensions(const hoNDArray<T>& x, hoNDArray<T>& r); // x : [... E1 E2], r: [... E2 E1]

    template<typename T> EXPORTCPUCOREMATH bool permuteE2To3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // x : [RO E1 CHA SLC E2 ...], r: [RO E1 E2 CHA SLC ...]
    template<typename T> EXPORTCPUCOREMATH bool permuteE2To5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // x : [RO E1 E2 CHA SLC ...], r: [RO E1 CHA SLC E2 ...]

    template<typename T> EXPORTCPUCOREMATH bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r); // x : [RO E1 E2 ...], r: [E1 E2 RO ...]
    template<typename T> EXPORTCPUCOREMATH bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r); // x : [RO E1 E2 CHA ...], r: [E1 E2 CHA RO ...]
    template<typename T> EXPORTCPUCOREMATH bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r); // x : [E1 E2 CHA RO ...], r: [RO E1 E2 CHA ...]

    template<typename T> EXPORTCPUCOREMATH bool permute3rdDimensionTo1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // x : [RO E1 E2 CHA ...], r: [E2 RO E1 CHA ...]

    template<typename T> EXPORTCPUCOREMATH bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r); // x : [RO E1 E2 srcCHA dstCHA ...], r: [E1 E2 srcCHA dstCHA RO ...]

    /// x : [RO E1 srcCHA], ker [RO E1 srcCHA dstCHA], buf is a buffer for computer, need to be pre-allocated [RO E1 srcCHA], y [RO E1 dstCHA]
    /// for the sake of speed, no check is made in this function
    template<typename T> EXPORTCPUCOREMATH bool imageDomainUnwrapping2D(const hoNDArray<T>& x, const hoNDArray<T>& ker, hoNDArray<T>& buf, hoNDArray<T>& y);

    /// x : [RO E1 srcCHA N], ker [RO E1 srcCHA dstCHA 1 or N], buf is a buffer for computer, need to be pre-allocated [RO E1 srcCHA], y [RO E1 dstCHA N]
    /// for the sake of speed, no check is made in this function
    template<typename T> EXPORTCPUCOREMATH bool imageDomainUnwrapping2DT(const hoNDArray<T>& x, const hoNDArray<T>& ker, hoNDArray<T>& buf, hoNDArray<T>& y);

#endif // USE_MKL
}
