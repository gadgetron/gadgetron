#pragma once

#include "hoCuNDArray.h"
#include "cudaDeviceManager.h"

namespace Gadgetron{

  template<class T> T dot( hoCuNDArray<T> *x, hoCuNDArray<T> *y, bool cc = true );
  
  template<class T> typename realType<T>::Type nrm2( hoCuNDArray<T> *x );
  
  template<class T> void axpy( T a, hoCuNDArray<T> *x, hoCuNDArray<T> *y );
  template<class T> void axpy( T a, hoCuNDArray< complext<T> > *x, hoCuNDArray< complext<T> > *y );
  
  /**
   * @brief Gets the index of the index of the element with minimum absolute
   * @param x Input data
   * @return index of absolute minimum values
   */
  template<class T> size_t amin( hoCuNDArray<T> *x);
  
  /**
   * @brief Gets the index of the index of the element with maximum absolute
   * @param x Input data
   * @return index of absolute maximum values
   * @details Note that this returns the C-style index and NOT the Fortran index.
   */
  template<class T> size_t amax( hoCuNDArray<T> *x );
  
  template<class T> typename realType<T>::Type asum( hoCuNDArray<T> *x );
}
