#pragma once

#include "cuNDArray.h"
#include "gpucore_export.h"

namespace Gadgetron{

  template<class T> EXPORTGPUCORE boost::shared_ptr<cuNDArray<T> > sum(cuNDArray<T> *data, unsigned int dim );
  
  template<class T> EXPORTGPUCORE T mean(cuNDArray<T>* data);
  
  template<class T> EXPORTGPUCORE T min(cuNDArray<T>* data);

  template<class T> EXPORTGPUCORE T max(cuNDArray<T>* data);
}
