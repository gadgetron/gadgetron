/** \file cuNDFFT.h
    \brief Wrapper of the CUFFT library for ndarrays of type Gadgetron::complext.
*/

#ifndef CUFFT_H
#define CUFFT_H
#pragma once

#include "cuNDArray.h"
#include "gpucore_export.h"

namespace Gadgetron{

  /** \class cuNDFFT
      \brief Wrapper of the CUFFT library for ndarrays of type complext.

      Wrapper of the CUFFT library for ndarrays of type complext<REAL>.
      The class' template type is a REAL, ie. float or double.
  */
  template<class T> class EXPORTGPUCORE cuNDFFT
  {
    
  public:
    cuNDFFT() {}
    virtual ~cuNDFFT() {}
    
    void fft ( cuNDArray<complext<T> > *input, std::vector<unsigned int> *dims_to_transform );
    void ifft( cuNDArray<complext<T> > *input, std::vector<unsigned int> *dims_to_transform, bool do_scale = true );
    
    void fft ( cuNDArray<complext<T> > *input, unsigned int dim_to_transform);
    void ifft( cuNDArray<complext<T> > *input, unsigned int dim_to_transform, bool do_scale = true );
    
    void fft ( cuNDArray<complext<T> > *input );
    void ifft( cuNDArray<complext<T> > *input, bool do_scale = true );
    
  protected:
    void fft_int( cuNDArray<complext<T> > *input, std::vector<unsigned int> *dims_to_transform, int direction, bool do_scale = true );
  };
}

#endif
