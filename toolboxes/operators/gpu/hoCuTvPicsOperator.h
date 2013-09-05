#pragma once


#include "hoCuTvOperator.h"
#include "tvPicsOperator.h"

namespace Gadgetron{

  template<class T, unsigned int D> class hoCuTvPicsOperator 
    : public tvPicsOperator< hoCuNDArray<T>, hoCuTvOperator<T,D>, typename realType<T>::Type >
  {
  public:
    hoCuTvPicsOperator() : tvPicsOperator< hoCuNDArray<T>, hoCuTvOperator<T,D>, typename realType<T>::Type >() {}
    virtual ~hoCuTvPicsOperator() {}
  };    
}
