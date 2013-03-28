#pragma once

#include "cuNDArray.h"
#include "vector_td.h"
#include "gpunfft_export.h"

namespace Gadgetron{

  // Crop
  template<class T, unsigned int D> EXPORTGPUNFFT
  void crop( typename uintd<D>::Type crop_offset, cuNDArray<T> *in, cuNDArray<T> *out );
  
  // Expand with zero filling
  template<class T, unsigned int D> EXPORTGPUNFFT
  void expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out );
  
  // Zero fill border (rectangular)
  template<class T, unsigned int D> EXPORTGPUNFFT
  void zero_fill_border( typename uintd<D>::Type matrix_size, cuNDArray<T> *image );
}
