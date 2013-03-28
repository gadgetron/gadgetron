#pragma once

#include "cuCgSolver.h"
#include "sbSolver.h"

#include "complext.h"

namespace Gadgetron{

  template <class T> class cuSbCgSolver : public sbSolver< cuNDArray<typename realType<T>::Type >, cuNDArray<T>, cuCgSolver<T> >
  {
  public:    
    cuSbCgSolver() : sbSolver<cuNDArray<typename realType<T>::Type >, cuNDArray<T>, cuCgSolver<T> >() {}    
    virtual ~cuSbCgSolver() {}
  };
}
