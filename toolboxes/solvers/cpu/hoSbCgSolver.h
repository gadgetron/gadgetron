#pragma once

#include "hoCgSolver.h"
#include "sbSolver.h"

#include "complext.h"

namespace Gadgetron{

  template <class T> class hoSbCgSolver : public sbSolver< hoNDArray<typename realType<T>::Type >, hoNDArray<T>, hoCgSolver<T> >
  {
  public:    
    hoSbCgSolver() : sbSolver<hoNDArray<typename realType<T>::Type >, hoNDArray<T>, hoCgSolver<T> >() {}    
    virtual ~hoSbCgSolver() {}
  };
}
