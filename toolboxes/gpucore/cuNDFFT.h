#ifndef CUNDFFT_H
#define CUNDFFT_H
#pragma once

#include "cuNDArray.h"
#include "gpucore_export.h"

/*
  Wrapper of the CUFFT library for ndarrays of complex type.
  ----------------------------------------------------------
  The class is instantiated for 
  - CUFFTs buildin complex types
  - Types complext (vector_td<REAL,2>)
*/

template<class T> class EXPORTGPUCORE cuNDFFT
{

 public:
  cuNDFFT() {}
  virtual ~cuNDFFT() {}

  void fft ( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform );
  void ifft( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform, bool do_scale = true );

  void fft ( cuNDArray<T> *input, unsigned int dim_to_transform);
  void ifft( cuNDArray<T> *input, unsigned int dim_to_transform, bool do_scale = true );

  void fft ( cuNDArray<T> *input );
  void ifft( cuNDArray<T> *input, bool do_scale = true );

 protected:
  void fft_int( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform, int direction, bool do_scale = true );
};

#endif
