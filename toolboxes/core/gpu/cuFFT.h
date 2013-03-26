/** \file cuFFT.h
    \brief Wrapper of the CUFFT library for ndarrays of complex type.
*/

#ifndef CUFFT_H
#define CUFFT_H
#pragma once

#include "cuNDArray.h"
#include "gpucore_export.h"

namespace Gadgetron{

  /** \class cuFFT
      \brief Wrapper of the CUFFT library for ndarrays of complex type.

      Wrapper of the CUFFT library for ndarrays of complex type.
      Supported (instantiated) types are
      - CUFFTs built-in complex types
      - complext<float> and complext<double>
  */
  template<class T> class EXPORTGPUCORE cuFFT
  {
    
  public:
    cuFFT() {}
    virtual ~cuFFT() {}
    
    void fft ( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform );
    void ifft( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform, bool do_scale = true );
    
    void fft ( cuNDArray<T> *input, unsigned int dim_to_transform);
    void ifft( cuNDArray<T> *input, unsigned int dim_to_transform, bool do_scale = true );
    
    void fft ( cuNDArray<T> *input );
    void ifft( cuNDArray<T> *input, bool do_scale = true );
    
  protected:
    void fft_int( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform, int direction, bool do_scale = true );
  };
}

#endif
