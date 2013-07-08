#pragma once

#include "cuNDArray.h"
#include "vector_td.h"
#include "gpucore_export.h"

namespace Gadgetron{

  template<class T> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
  shift_dim( cuNDArray<T> *in, int shift );

  template<class T> EXPORTGPUCORE void
  shift_dim( cuNDArray<T> *in, cuNDArray<T> *out, int shift );
  
  template<class T> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
  permute( cuNDArray<T> *in, std::vector<unsigned int> *dim_order, int shift_mode = 0 );
  
  template<class T> EXPORTGPUCORE void
  permute( cuNDArray<T> *in, cuNDArray<T> *out, std::vector<unsigned int> *dim_order, int shift_mode = 0 );

  template<class T, unsigned int D> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
  crop( typename uintd<D>::Type crop_offset, typename uintd<D>::Type crop_size, cuNDArray<T> *in );

  template<class T, unsigned int D> EXPORTGPUCORE
  void crop( typename uintd<D>::Type crop_offset, cuNDArray<T> *in, cuNDArray<T> *out );
  
  template<class T, unsigned int D> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
  pad( typename uintd<D>::Type size, cuNDArray<T> *in, T val = T(0) );

  template<class T, unsigned int D> EXPORTGPUCORE
  void pad( cuNDArray<T> *in, cuNDArray<T> *out, T val = T(0) );
  
  template<class T, unsigned int D> EXPORTGPUCORE
  void fill_border( typename uintd<D>::Type matrix_size, cuNDArray<T> *image, T val = T(0) );

  // Expand array to new dimension
  template<class T> EXPORTGPUCORE boost::shared_ptr<cuNDArray<T> > 
  expand(cuNDArray<T> *data, unsigned int added_dim_size );
  
  // Sum over dimension
  template<class T> EXPORTGPUCORE boost::shared_ptr<cuNDArray<T> > 
  sum(cuNDArray<T> *data, unsigned int dim );

  template<class T> EXPORTGPUCORE T mean(cuNDArray<T>* data);
}
