#pragma once

#include "sbSolver.h"
#include "cuCgSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"

#include "complext.h"

namespace Gadgetron{
template <class T> class cuSbCgSolver
  : public sbSolver<cuNDArray<typename realType<T>::type >, cuNDArray<T>, cuCgSolver<T> >
{
public:
  
  cuSbCgSolver() : sbSolver<cuNDArray<typename realType<T>::type >, cuNDArray<T>, cuCgSolver<T> >() {
  }
  
  virtual ~cuSbCgSolver() {}
};
}
