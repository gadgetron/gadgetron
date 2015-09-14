#pragma once

#include "hoCuNDArray_math.h"
#include "encodingOperatorContainer.h"

namespace Gadgetron{
  
  template<class T> class hoCuEncodingOperatorContainer 
    : public encodingOperatorContainer< hoCuNDArray<T> >
  {
  public:
    hoCuEncodingOperatorContainer() : encodingOperatorContainer< hoCuNDArray<T> >() {}
    virtual ~hoCuEncodingOperatorContainer() {}
    
  }; 
}
