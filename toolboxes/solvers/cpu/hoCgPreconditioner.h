#pragma once

#include "hoNDArray_math.h"
#include "cgPreconditioner.h"

namespace Gadgetron{

  template<class T> class hoCgPreconditioner : public cgPreconditioner< hoNDArray<T> >
  {
  public:    
    hoCgPreconditioner() : cgPreconditioner< hoNDArray<T> >() {}
    virtual ~hoCgPreconditioner() {}
  };
}
