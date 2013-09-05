#pragma once

#include "hoNDArray_math.h"
#include "hoTvOperator.h"
#include "tvPicsOperator.h"

namespace Gadgetron{

  template<class T, unsigned int D> class hoTvPicsOperator 
    : public tvPicsOperator< hoNDArray<T>, hoTvOperator<T,D>, typename realType<T>::Type >
  {
  public:
    hoTvPicsOperator() : tvPicsOperator< hoNDArray<T>, hoTvOperator<T,D>, typename realType<T>::Type >() {}
    virtual ~hoTvPicsOperator() {}
  };    
}
