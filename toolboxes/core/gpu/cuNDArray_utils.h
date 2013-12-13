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
  permute( cuNDArray<T> *in, std::vector<size_t> *dim_order, int shift_mode = 0 );
  
  template<class T> EXPORTGPUCORE void
  permute( cuNDArray<T> *in, cuNDArray<T> *out, std::vector<size_t> *dim_order, int shift_mode = 0 );

  template<class T, unsigned int D> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
  crop( typename uint64d<D>::Type crop_offset, typename uint64d<D>::Type crop_size, cuNDArray<T> *in );

  template<class T, unsigned int D> EXPORTGPUCORE
  void crop( typename uint64d<D>::Type crop_offset, cuNDArray<T> *in, cuNDArray<T> *out );
  
  template<class T, unsigned int D> EXPORTGPUCORE boost::shared_ptr< cuNDArray<T> >
  pad( typename uint64d<D>::Type size, cuNDArray<T> *in, T val = T(0) );

  template<class T, unsigned int D> EXPORTGPUCORE
  void pad( cuNDArray<T> *in, cuNDArray<T> *out, T val = T(0) );
  
  template<class T, unsigned int D> EXPORTGPUCORE
  void fill_border( typename uint64d<D>::Type matrix_size, cuNDArray<T> *image, T val = T(0) );

  /***
   * @brief Fills the image with a given value outside a radius from the center
   * @param radius
   * @param in_out
   * @param val
   */
  template<class T, unsigned int D>
  void fill_border( typename realType<T>::Type radius, cuNDArray<T> *in_out, T val= T(0) );

  // Expand array to new dimension
  template<class T> EXPORTGPUCORE boost::shared_ptr<cuNDArray<T> > 
  expand(cuNDArray<T> *data, size_t added_dim_size );
  
  template<class T, unsigned int D> EXPORTGPUCORE 
  boost::shared_ptr< cuNDArray<T> > upsample( cuNDArray<T>* in );

  template<class T, unsigned int D> EXPORTGPUCORE
  void upsample( cuNDArray<T> *in, cuNDArray<T> *out );

  template<class T, unsigned int D> EXPORTGPUCORE 
  boost::shared_ptr< cuNDArray<T> > downsample( cuNDArray<T>* in );

  template<class T, unsigned int D> EXPORTGPUCORE
  void downsample( cuNDArray<T> *in, cuNDArray<T> *out );
}
