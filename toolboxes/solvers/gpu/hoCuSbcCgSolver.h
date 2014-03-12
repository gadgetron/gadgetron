#pragma once

#include "hoCuCgSolver.h"
#include "sbcSolver.h"

#include "complext.h"

namespace Gadgetron{

  template <class T> class hoCuSbcCgSolver : public sbcSolver< hoCuNDArray<typename realType<T>::Type >, hoCuNDArray<T>, hoCuCgSolver<T> >
  {
  public:
    hoCuSbcCgSolver() : sbcSolver<hoCuNDArray<typename realType<T>::Type >, hoCuNDArray<T>, hoCuCgSolver<T> >() {}
    virtual ~hoCuSbcCgSolver() {}
  };
}
