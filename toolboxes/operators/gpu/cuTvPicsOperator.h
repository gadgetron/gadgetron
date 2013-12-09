#pragma once

#include "cuNDArray.h"
#include "cuTvOperator.h"
#include "tvPicsOperator.h"

namespace Gadgetron{

  template<class T, unsigned long long D> class cuTvPicsOperator 
    : public tvPicsOperator< cuNDArray<T>, cuTvOperator<T,D>, typename realType<T>::Type >
  {
  public:
    cuTvPicsOperator() : tvPicsOperator< cuNDArray<T>, cuTvOperator<T,D>, typename realType<T>::Type >() {}
    virtual ~cuTvPicsOperator() {}
  };    
}
