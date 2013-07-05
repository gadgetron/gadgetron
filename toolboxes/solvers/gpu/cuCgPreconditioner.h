#pragma once

#include "cuNDArray_operators.h"
#include "cgPreconditioner.h"

namespace Gadgetron{

  template<class T> class cuCgPreconditioner : public cgPreconditioner< cuNDArray<T> >
  {
  public:    
    cuCgPreconditioner() : cgPreconditioner< cuNDArray<T> >() {}
    virtual ~cuCgPreconditioner() {}
  };
}
