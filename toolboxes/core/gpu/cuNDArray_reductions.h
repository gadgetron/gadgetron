#pragma once

#include "cuNDArray.h"

namespace Gadgetron{

  template<class T> boost::shared_ptr<cuNDArray<T> > sum(cuNDArray<T> *data, unsigned int dim );
  
  template<class T> T mean(cuNDArray<T>* data);
  
  template<class T> T min(cuNDArray<T>* data);

  template<class T> T max(cuNDArray<T>* data);
}
