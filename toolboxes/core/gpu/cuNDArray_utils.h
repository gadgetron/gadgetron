#pragma once

#include "cuNDArray.h"
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

  // Expand array to new dimension
  template<class T> EXPORTGPUCORE boost::shared_ptr<cuNDArray<T> > expand(cuNDArray<T> *data, unsigned int added_dim_size );
  
  // Sum over dimension
  template<class T> EXPORTGPUCORE boost::shared_ptr<cuNDArray<T> > sum(cuNDArray<T> *data, unsigned int dim );
}
