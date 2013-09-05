#pragma once

#include "cuNDArray_math.h"
#include "cuTvOperator.h"
#include "tvPicsOperator.h"

namespace Gadgetron{

  template<class T, unsigned int D> class cuTvPicsOperator 
    : public tvPicsOperator< cuNDArray<T>, cuTvOperator<T,D>, typename realType<T>::Type >
  {
  public:
    cuTvPicsOperator() : tvPicsOperator< cuNDArray<T>, cuTvOperator<T,D>, typename realType<T>::Type >() {}
    virtual ~cuTvPicsOperator() {}
  };    
}
