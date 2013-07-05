#pragma once

#include "cuCgSolver.h"
#include "sbcSolver.h"

namespace Gadgetron{
  
  template <class T> class cuSbcCgSolver : public sbcSolver< cuNDArray<typename realType<T>::Type >, cuNDArray<T>, cuCgSolver<T> >
  {
  public:    
    cuSbcCgSolver() : sbcSolver<cuNDArray<typename realType<T>::Type >, cuNDArray<T>, cuCgSolver<T> >() {}
    virtual ~cuSbcCgSolver() {}    
  };
}
