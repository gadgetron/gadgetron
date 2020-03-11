/** \file   hoNDArray_elemwise.h
    \brief  Element-wise math operations on the hoNDArray class.

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

    3) Many functions are also reimplemented if the Intel MKL is avaiable to speedup the computation.
 */

#pragma once

#include "hoNDArray.h"
#include "cpucore_math_export.h"
#include "cpp_blas.h"

#include <complex>

namespace Gadgetron{

  //
  // Math return types
  //
  template <class T, class I> struct mathReturnType {};

  template <class T> struct mathReturnType<T,T> {typedef T type;};

  template <class T> struct mathReturnType<complext<T>,T> {typedef complext<T> type;};
  template <class T> struct mathReturnType<T,complext<T> > {typedef complext<T> type;};
  template <class T> struct mathReturnType<complext<T>,complext<T> > {typedef complext<T> type;};

  template <class T> struct mathReturnType<std::complex<T>,T> {typedef std::complex<T> type;};
  template <class T> struct mathReturnType<T,std::complex<T> > {typedef std::complex<T> type;};
  template <class T> struct mathReturnType<std::complex<T>,std::complex<T> > {typedef std::complex<T> type;};

  template <class T, class S> struct mathReturnType<T, complext<S> > {typedef complext<typename mathReturnType<T,S>::type> type;};
  template <class T, class S> struct mathReturnType<complext<T>, S> {typedef complext<typename mathReturnType<T,S>::type> type;};
  template <class T, class S> struct mathReturnType<complext<T>, complext<S> > {typedef complext<typename mathReturnType<T,S>::type> type;};

  template <class T, class S> struct mathReturnType<T, std::complex<S> > {typedef std::complex<typename mathReturnType<T,S>::type> type;};
  template <class T, class S> struct mathReturnType<std::complex<T>, S> {typedef std::complex<typename mathReturnType<T,S>::type> type;};
  template <class T, class S> struct mathReturnType<std::complex<T>, std::complex<S> > {typedef std::complex<typename mathReturnType<T,S>::type> type;};

  template<> struct mathReturnType<unsigned int, int> {typedef int type;};
  template<> struct mathReturnType<int, unsigned int> {typedef int type;};
  template<> struct mathReturnType<int, bool> {typedef int type;};
  template<> struct mathReturnType<bool,int> {typedef int type;};
  template<> struct mathReturnType<unsigned int, bool> {typedef int type;};
  template<> struct mathReturnType<bool,unsigned int> {typedef int type;};
  template<> struct mathReturnType<float, unsigned int> {typedef float type;};
  template<> struct mathReturnType<unsigned int, float> {typedef float type;};
  template<> struct mathReturnType<float, int> {typedef float type;};
  template<> struct mathReturnType<int, float> {typedef float type;};
  template<> struct mathReturnType<float, bool> {typedef float type;};
  template<> struct mathReturnType<bool, float> {typedef float type;};
  template<> struct mathReturnType<double, unsigned int> {typedef double type;};
  template<> struct mathReturnType<unsigned int, double> {typedef double type;};
  template<> struct mathReturnType<double, int> {typedef double type;};
  template<> struct mathReturnType<int, double> {typedef double type;};
  template<> struct mathReturnType<double, bool> {typedef double type;};
  template<> struct mathReturnType<bool, double> {typedef double type;};
  template<> struct mathReturnType<double, float> {typedef double type;};
  template<> struct mathReturnType<float,double> {typedef double type;};

  // Utility to verify array dimensions for simple broadcasting.
  // It "replaces" NDArray::dimensions_equal() to support batch mode.
  // There is an identical function for all array instances (currently hoNDArray, cuNDArray, hoCuNDArray)
  // !!! Remember to fix any bugs in all versions !!!
  //
  template<class T,class S> bool compatible_dimensions( const hoNDArray<T> &x, const hoNDArray<S> &y )
  {
      return ((x.get_number_of_elements()%y.get_number_of_elements())==0);
  }

  // Utility to verify if array dimensions are compatible for a binary operation
  // that supports simple broadcasting, i.e. for f(x,y,r), there are three cases:
  // 1) nr = nx = ny
  // 2) nr = nx > ny, and nx is divisible by ny
  // 3) nr = ny > nx, and ny is divisible by nx
  //
  template<class T,class S, class U> bool compatible_dimensions( const hoNDArray<T> &x, const hoNDArray<S> &y,
              const hoNDArray<U> &r )
  {
      size_t nx = x.get_number_of_elements();
      size_t ny = y.get_number_of_elements();
      size_t nr = r.get_number_of_elements();
      if (nx == ny) {
          return (nx==nr);
      }
      if ((nx%ny)==0) {
          return (nx==nr);
      }
      if ((ny%nx)==0) {
          return (ny==nr);
      }
      return false;
  }


/**
* @brief add two vectors of values, r = x + y
  support in-place computation, e.g. x==r or y==r
  support simple broadcasting
*/
template <class T, class S>
void add(const hoNDArray<T>& x, const hoNDArray<S>& y, hoNDArray<typename mathReturnType<T,S>::type >& r);

// Pointer version calls the reference version
template <typename T, class S>
void add(const hoNDArray<T>* x, const hoNDArray<S>* y, hoNDArray<typename mathReturnType<T,S>::type >* r)
{
  add(*x, *y, *r);
}


/**
* @brief subtract two vectors of values, r = x - y
  support in-place computation, e.g. x==r
  support simple broadcasting
*/
template <class T, class S>
void subtract(const hoNDArray<T>& x, const hoNDArray<S>& y, hoNDArray<typename mathReturnType<T,S>::type >& r);

// Pointer version calls the reference version
template <typename T, class S>
void subtract(const hoNDArray<T>* x, const hoNDArray<S>* y, hoNDArray<typename mathReturnType<T,S>::type >* r)
{
  subtract(*x, *y, *r);
}


/**
* @brief multiply two vectors of values, r = x * y
  support in-place computation, e.g. x==r or y==r
  support simple broadcasting
*/
template <class T, class S>
void multiply(const hoNDArray<T>& x, const hoNDArray<S>& y, hoNDArray<typename mathReturnType<T,S>::type >& r);

// Pointer version calls the reference version
template <class T, class S>
void multiply(const hoNDArray<T>* x, const hoNDArray<S>* y, hoNDArray<typename mathReturnType<T,S>::type >* r)
{
  multiply(*x, *y, *r);
}


/**
* @brief divide two vectors of values, r = x / y
  support in-place computation, e.g. x==r
  support simple broadcasting
  no check for y==0
*/
template <class T, class S>
void divide(const hoNDArray<T>& x, const hoNDArray<S>& y, hoNDArray<typename mathReturnType<T,S>::type >& r);

// Pointer version calls the reference version
template <class T, class S>
void divide(const hoNDArray<T>* x, const hoNDArray<S>* y, hoNDArray<typename mathReturnType<T,S>::type >* r)
{
  divide(*x, *y, *r);
}


/**
* @brief r = x * conj(y)
  support in-place computation, e.g. x==r
  support simple broadcasting
*/
template <class T, class S>
void multiplyConj(const hoNDArray<T>& x, const hoNDArray<S>& y, hoNDArray<typename mathReturnType<T,S>::type >& r);

// Pointer version calls the reference version
template <class T, class S>
void multiplyConj(const hoNDArray<T>* x, const hoNDArray<S>* y, hoNDArray<typename mathReturnType<T,S>::type >* r)
{
  multiplyConj(*x, *y, *r);
}


/**
* @brief r = conj(x)
*/
template <typename T>
void conjugate(const hoNDArray<T>& x, hoNDArray<T>& r);

/**
* @brief if abs(x) is smaller than epsilon for its numeric type
add epsilon to this x
*/
template <typename T>
void addEpsilon(hoNDArray<T>& x);

/**
* @brief r = angle(x)
*/
template <typename T>
void argument(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r);
template<class T>
hoNDArray<realType_t<T>> argument(const hoNDArray<T>& x);

/**
* @brief r = 1/x
*/
template <typename T>
void inv(const hoNDArray<T>& x, hoNDArray<T>& r);

/**
 * @brief Calculates the element-wise absolute values (l2 norm) of the array entries
 * @param[in] x Input array.
 * @return A new array containing the element-wise absolute values of the input.
 */
template<class T> boost::shared_ptr< hoNDArray<typename realType<T>::Type> > abs( hoNDArray<T> *x );
template<class T> hoNDArray<typename realType<T>::Type> abs(const  hoNDArray<T>& x );
template <typename T, typename R> void abs(const hoNDArray<T>& x, hoNDArray<R>& r);

/**
 * @brief Calculates the element-wise absolute values (l2 norm) of the array entries (in place).
 * @param[in,out] x Input and output array.
 */
template<class T> void abs_inplace( hoNDArray<T> *x );

/**
 * @brief Calculates the element-wise squared absolute values of the array entries
 * @param[in] x Input array.
 * @return A new array containing the element-wise absolute values of the input.
 */
template<class T> boost::shared_ptr< hoNDArray<typename realType<T>::Type> > abs_square( hoNDArray<T> *x );

/**
 * @brief Calculates the element-wise sqrt of the array entries.
 * @param[in] x Input array.
 * @return A new array containing the element-wise sqrt of the input.
 */
template<class T> boost::shared_ptr< hoNDArray<T> > sqrt( hoNDArray<T> *x );

template <typename T> void sqrt(const hoNDArray<T>& x, hoNDArray<T>& r);

/**
 * @brief Calculates the element-wise sqrt of the array entries (in place).
 * @param[in,out] x Input and output array.
 */
template<class T> void sqrt_inplace( hoNDArray<T> *x );

/**
 * @brief Calculates the element-wise square of the array entries.
 * @param[in] x Input array.
 * @return A new array containing the element-wise square of the input.
 *
 * For real numbers this functions is equivalent to square.
 * For complex arrays abs_square() and square() differ however.
 */
template<class T> boost::shared_ptr< hoNDArray<T> > square( hoNDArray<T> *x );

/**
 * @brief Calculates the element-wise square of the array entries (in place).
 * @param[in,out] x Input and output array.
 */
template<class T> void square_inplace( hoNDArray<T> *x );

/**
 * @brief Calculates the element-wise reciprocal of the array entries.
 * @param[in] x Input array.
 * @return A new array containing the element-wise reciprocal of the input.
 */
template<class T> boost::shared_ptr< hoNDArray<T> > reciprocal( hoNDArray<T> *x );

/**
 * @brief Calculates the element-wise reciprocal of the array entries (in place).
 * @param[in,out] x Input and output array.
 */
template<class T> void reciprocal_inplace( hoNDArray<T> *x );

/**
 * @brief Calculates the element-wise reciprocal sqrt of the array entries.
 * @param[in] x Input array.
 * @return A new array containing the element-wise reciprocal sqrt of the input.
 */
template<class T> boost::shared_ptr< hoNDArray<T> > reciprocal_sqrt( hoNDArray<T> *x );

/**
 * @brief Calculates the element-wise reciprocal sqrt of the array entries (in place).
 * @param[in,out] x Input and output array.
 */
template<class T> void reciprocal_sqrt_inplace( hoNDArray<T> *x );

/**
 * @brief Calculates the elementwise signum function on the array.
 * @param[in] x Input array.
 * @return A new array containing the element-wise sgn of the input.
 */
template<class T> boost::shared_ptr< hoNDArray<T> > sgn( hoNDArray<T> *x );

/**
 * @brief Calculates the elementwise signum function on the array (in place).
 * @param[in,out] x Input and output array.
 */
template<class T> void sgn_inplace( hoNDArray<T> *x );

/**
 * @brief Extract the real component from a complex array.
 * @param[in] x Input array.
 * @return A new array of the real component of the complex array.
 */
template<class T> boost::shared_ptr< hoNDArray<typename realType<T>::Type> > real(const  hoNDArray<T> *x );
    template<class T>  hoNDArray<realType_t<T>> real(const  hoNDArray<T>&x );

/**
 * @brief Extract the imaginary component from a complex array.
 * @param[in] x Input array.
 * @return A new array of the imaginary component of the complex array.
 */
template<class T> boost::shared_ptr< hoNDArray<typename realType<T>::Type> > imag( const hoNDArray<T> *x );
template<class T>  hoNDArray<realType_t<T>> imag(const  hoNDArray<T>&x );
/**
 * @brief Create a new array of the complex conjugate of the input array. For real arrays a copy of the input array is return.
 * @param[in] x Input array.
 * @return A new array of the complex conjugate of the input array.
 */
template<class T> boost::shared_ptr< hoNDArray<T> > conj( const hoNDArray<T> *x );

/**
 * @brief Construct a complex array from a real array.
 * @param[in] x Input array.
 * @return A new complex array containing the input array in the real component and zeros in the imaginary component.
 */
template<class T> boost::shared_ptr< hoNDArray<T> >
real_to_complex( const hoNDArray<typename realType<T>::Type> *x );

template<class T> boost::shared_ptr< hoNDArray<T> >
real_imag_to_complex( hoNDArray<typename realType<T>::Type> *real, hoNDArray<typename realType<T>::Type>* imag);

/**
* @brief real and imag to complex
*/
template<class T>
void real_imag_to_complex(const hoNDArray<typename realType<T>::Type>& real, const hoNDArray<typename realType<T>::Type>& imag, hoNDArray<T>& cplx);

/**
* @brief complex to real and imag
*/
template<class T>
void complex_to_real_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real, hoNDArray<typename realType<T>::Type>& imag);

template<class T>
void complex_to_real_imag(const hoNDArray<T>& cplx, hoNDArray<T>& real, hoNDArray<T>& imag);

/**
* @brief get the real part of complex
*/
template<class T>
void complex_to_real(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real);

template<class T>
void complex_to_real(const hoNDArray<T>& cplx, hoNDArray<T>& real);

template<class T> 
void complex_to_real(hoNDArray<T>& cplx);

/**
* @brief get the imag part of complex
*/
template<class T>
void complex_to_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& imag);

template<class T>
void complex_to_imag(const hoNDArray<T>& cplx, hoNDArray<T>& imag);

template<class T>
void complex_to_imag(hoNDArray<T>& cplx);

/**
* @brief get complex array whose real part is the input and imag part is zero
*/
template<class T>
void real_to_complex(const hoNDArray<typename realType<T>::Type>& real, hoNDArray<T>& cplx);

/**
 * @brief Clears the array to all zeros ( in place). Faster than fill.
 * @param[in,out] x Input and output array.
 */
template<class T> void clear( hoNDArray<T>* x )
{
    if ( x->get_number_of_elements() > 0 )
    {
        memset( x->get_data_ptr(), 0, x->get_number_of_elements()*sizeof(T));
    }
}

template<class T> void clear( hoNDArray<T>& x )
{
    if ( x.get_number_of_elements() > 0 )
    {
        memset( x.get_data_ptr(), 0, x.get_number_of_elements() * sizeof(T));
    }
}

/**
 * @brief Fills the array with a user provided constant value (in place).
 * @param[in,out] x Input and output array.
 * @param[in] val Fill value.
 */
template <typename T> void fill( hoNDArray<T>* x, T val);
template <typename T> void fill( hoNDArray<T>& x, T val );

/**
 * @brief Clamps all values in the array to the minimum and maximum values specified (in place).
 * @param[in,out] x Input and output array.
 * @param[in] min minimum value.
 * @param[in] max maximum value.
 * @param[in] min_val value to which everything below the minimum will be set
 * @param[in] max_val value to which everything above the maximum will be set
 */
template<class T> void clamp( hoNDArray<T> *x, typename realType<T>::Type min, typename realType<T>::Type max, T min_val, T max_val );

/**
 * @brief Clamps all values in the array to the minimum and maximum values specified (in place).
 * @param[in,out] x Input and output array.
 * @param[in] min minimum value.
 * @param[in] max maximum value.
 */
template<class T> void clamp( hoNDArray<T> *x, typename realType<T>::Type min, typename realType<T>::Type max );

/**
 * @brief Clamps all values in the array to a minimum value allowed (in place).
 * @param[in,out] x Input and output array.
 * @param[in] min Minimum value.
 */
template<class T> void clamp_min( hoNDArray<T> *x, typename realType<T>::Type min );

/**
 * @brief Clamps all values in the array to a maximum value allowed (in place).
 * @param[in,out] x Input and output array.
 * @param[in] max Maximum value.
 */
template<class T> void clamp_max( hoNDArray<T> *x, typename realType<T>::Type max );

/**
 * @brief In place normalization (scaling) to a new maximum absolute array value val.
 * @param[in,out] x Input and output array.
 * @param[in] val New maximum absolute array value (according to the l2-norm)
 */
template<class T> void normalize( hoNDArray<T> *x, typename realType<T>::Type val = typename realType<T>::Type(1) );

/**
 * @brief Shrinkage (soft thresholding), i.e. shrink(x,gamma) = x/abs(x)*max(abs(x)-gamma,0).
 * @param[out] out Output array. Can be 0x0 in which case an in place transform is performed.
 * @param[in,out] x Input array (and output array if out == 0x0).
 * @param[in] gamma Shrinkage control parameter
 */
template<class T> void shrink1( hoNDArray<T> *x, typename realType<T>::Type gamma, hoNDArray<T> *out = 0x0 );

/**
 * @brief In place p-shrinkage (soft thresholding), i.e. pshrink(x,gamma,p) = x/abs(x)*max(abs(x)-gamma*abs(x)^(p-1),0).
 * @param[out] out Output array. Can be 0x0 in which case an in place transform is performed.
 * @param[in,out] x Input array (and output array if out == 0x0).
 * @param[in] gamma Shrinkage control parameter
 * @param[in] p p value of the shrinkage. Should be less than 1 and more than 0.
 */
template<class T> void pshrink( hoNDArray<T> *x, typename realType<T>::Type gamma,typename realType<T>::Type p, hoNDArray<T> *out = 0x0 );

/**
 * @brief Shrinkage (soft thresholding, multi-dimensional), i.e. shrink(x,gamma,s) = x/s*max(s-gamma,0).
 * @param[out] out Output array. Can be 0x0 in which case an in place transform is performed.
 * @param[in,out] x Input array (and output array if out == 0x0).
 * @param[in] s Input array, normalization.
 * @param[in] gamma Shrinkage control parameter
 */
template<class T> void shrinkd ( hoNDArray<T> *x, hoNDArray<typename realType<T>::Type> *s, typename realType<T>::Type gamma, hoNDArray<T> *out = 0x0 );

/**
 * @brief In place p-shrinkage (soft thresholding, multi-dimensional), i.e. pshrink(x,s,gamma,p) = x/s*max(s-gamma*s^(p-1),0).
 * @param[out] out Output array. Can be 0x0 in which case an in place transform is performed.
 * @param[in,out] x Input array (and output array if out == 0x0).
 * @param[in] gamma Shrinkage control parameter
 * @param[in] p p value of the shrinkage. Should be less than 1 and more than 0.
 */
template<class T> void pshrinkd ( hoNDArray<T> *x, hoNDArray<typename realType<T>::Type> *s, typename realType<T>::Type gamma,typename realType<T>::Type p, hoNDArray<T> *out = 0x0 );

/**
 * @brief Implementation of element-wise operator+= on two hoNDArrays.
 * @param[in,out] x Input and output array.
 * @param[in] y Input array.

 * Let y be an n-dimensional array.
 * Then the sizes of the first n array dimensions must match between x and y.
 * If x contains further dimensions the operator is batched across those dimensions.
 */
template<class T, class S> hoNDArray<T>& operator+= (hoNDArray<T> &x, const hoNDArray<S> &y);

/**
 * @brief Implementation of element-wise operator+= on a hoNDArray with a scalar value.
 * @param[in,out] x Input and output array.
 * @param[in] y Input scalar.
 */
template<class T, class R> hoNDArray<T>& operator+= (hoNDArray<T> &x, const R &y);


/**
 * @brief Implementation of element-wise operator-= on two hoNDArrays.
 * @param[in,out] x Input and output array.
 * @param[in] y Input array.

 * Let y be an n-dimensional array.
 * Then the sizes of the first n array dimensions must match between x and y.
 * If x contains further dimensions the operator is batched across those dimensions.
 */
template<class T, class S> hoNDArray<T>& operator-= (hoNDArray<T> &x, const hoNDArray<S> &y);


/**
 * @brief Implementation of element-wise operator-= on a hoNDArray with a scalar value.
 * @param[in,out] x Input and output array.
 * @param[in] y Input scalar.
 */
template<class T, class S> hoNDArray<T>& operator-= (hoNDArray<T> &x, const S &y);


/**
 * @brief Implementation of element-wise operator*= on two hoNDArrays.
 * @param[in,out] x Input and output array.
 * @param[in] y Input array.

 * Let y be an n-dimensional array.
 * Then the sizes of the first n array dimensions must match between x and y.
 * If x contains further dimensions the operator is batched across those dimensions.
 */
template<class T, class S> hoNDArray<T>& operator*= (hoNDArray<T> &x, const hoNDArray<S> &y);

/**
 * @brief Implementation of element-wise operator*= on a hoNDArray with a scalar value.
 * @param[in,out] x Input and output array.
 * @param[in] y Input scalar.
 */
template<class T, class S> hoNDArray<T>& operator*= (hoNDArray<T> &x, const S &y);

/**
 * @brief Implementation of element-wise operator/= on two hoNDArrays.
 * @param[in,out] x Input and output array.
 * @param[in] y Input array.

 * Let y be an n-dimensional array.
 * Then the sizes of the first n array dimensions must match between x and y.
 * If x contains further dimensions the operator is batched across those dimensions.
 */
template<class T, class S> hoNDArray<T>& operator/= (hoNDArray<T> &x, const hoNDArray<S> &y);

/**
 * @brief Implementation of element-wise operator/= on a hoNDArray with a scalar value.
 * @param[in,out] x Input and output array.
 * @param[in] y Input scalar.
 */
template<class T, class S> hoNDArray<T>& operator/= (hoNDArray<T> &x, const S &y);




/**
 * @brief Calculates y = a*x+y in which x and y are considered as vectors
 * @param[in] a Scalar value
 * @param[in] x Array
 * @param[in,out] y Array
 */
template<class T> void axpy(T a, const hoNDArray<T> *x, hoNDArray<T> *y ){ BLAS::axpy(x->get_number_of_elements(),a,x->get_data_ptr(),1,y->get_data_ptr(),1);}

/**
* @brief compute r = a*x + y
*/
template <typename T> void axpy(T a, const hoNDArray<T>& x, hoNDArray<T>& y){ axpy(a,&x,&y);}

/**
* @brief compute x *= a
*/
template <typename R, typename T> void scal(R a, hoNDArray<T>& x) {BLAS::scal(x.get_number_of_elements(),a,x.get_data_ptr(),1);}

/**
* @brief 2D convolution
            x: input data, y: convolution kernel, z: output; each 2D slice is convolved
*/
template <typename T>
void conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z);

/**
* @brief 3D convolution
            x: input data, y: convolution kernel, z: output; each 3D volume is convolved
*/
template <typename T>
void conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z);

/**
* @brief sum over a specific dimension
            x: input array, y: output array, dim: dimension to perform sum
            resulting y.get_size(d) == 1
*/
template <typename T>
void sum_over_dimension(const hoNDArray<T>& x, hoNDArray<T>& y, size_t dim);



/**
 * @brief Implementation of element-wise operator&= on two hoNDArrays.
 * @param[in,out] x Input and output array.
 * @param[in] y Input array.

 * Let y be an n-dimensional array.
 * Then the sizes of the first n array dimensions must match between x and y.
 * If x contains further dimensions the operator is batched across those dimensions.
 */
 hoNDArray<bool>& operator&= (hoNDArray<bool> &x, const hoNDArray<bool> &y);

/**
 * @brief Implementation of element-wise operator&= on two hoNDArrays.
 * @param[in,out] x Input and output array.
 * @param[in] y Input array.

 * Let y be an n-dimensional array.
 * Then the sizes of the first n array dimensions must match between x and y.
 * If x contains further dimensions the operator is batched across those dimensions.
 */
 hoNDArray<bool>& operator|= (hoNDArray<bool> &x, const hoNDArray<bool> &y);

 /**
  * Function transforming an input array into an output array
  * @tparam T Element type of input array
  * @tparam S Element type of output array
  * @tparam F Function type
  * @param input Input array
  * @param output Output array. Must have the same number of elements as the input array.
  * @param fun Function taking a value of type T and returning a value assignable to type S.
  */
 template <class T, class S, class F> void transform(const hoNDArray<T>& input, hoNDArray<S>& output, F&& fun);


 /**
  * Returns a new array containing the input array transformed by the provided function
  * @tparam T Element type of input array
  * @tparam F Function type
  * @tparam S Element type of output array.
  * @param input  Input array to be transformed
  * @param fun Function taking a value of type T and returning a value assignable to type S.
  * @return The resulting array
  */
 template <class T, class F, class S=std::invoke_result_t<F&&,T&&>> hoNDArray<S> transform(const hoNDArray<T>& input, F&& fun);

}


#include "hoNDArray_elemwise.hpp"
