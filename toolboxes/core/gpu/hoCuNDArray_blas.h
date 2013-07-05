#pragma once
#include "hoCuNDArray.h"
#include "cudaDeviceManager.h"
#include "gpucore_export.h"

namespace Gadgetron{

  template<class T> EXPORTGPUCORE T dot( hoCuNDArray<T> *x, hoCuNDArray<T> *y, bool cc = true );
  
  template<class T> EXPORTGPUCORE typename realType<T>::Type nrm2( hoCuNDArray<T> *x );
  
  template<class T> EXPORTGPUCORE void axpy( T a, hoCuNDArray<T> *x, hoCuNDArray<T> *y );
  template<class T> EXPORTGPUCORE void axpy( T a, hoCuNDArray< complext<T> > *x, hoCuNDArray< complext<T> > *y );
  
  /**
   * @brief Gets the index of the index of the element with minimum absolute
   * @param x Input data
   * @return index of absolute minimum values
   */
  template<class T> EXPORTGPUCORE int amin( hoCuNDArray<T> *x);
  
  /**
   * @brief Gets the index of the index of the element with maximum absolute
   * @param x Input data
   * @return index of absolute maximum values
   * @details Note that this returns the C-style index and NOT the Fortran index.
   */
  template<class T> EXPORTGPUCORE int amax( hoCuNDArray<T> *x );
  
  template<class T> EXPORTGPUCORE typename realType<T>::Type asum( hoCuNDArray<T> *x );
}
