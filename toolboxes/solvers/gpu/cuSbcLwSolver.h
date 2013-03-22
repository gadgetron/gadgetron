#pragma once

#include "sbcSolver.h"
#include "cuLwSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "encodingOperatorContainer.h"

namespace Gadgetron{
template <class T> class cuSbcLwSolver
  : public sbcSolver< cuNDArray<typename realType<T>::type>, cuNDArray<T>, cuLwSolver<T> >
{
public:
  
  cuSbcLwSolver() : sbcSolver<cuNDArray<typename realType<T>::type>, cuNDArray<T>, cuLwSolver<T> >() {
    set_device(-1); 
  }

  virtual ~cuSbcLwSolver() {}

#include "cuSbSolver_macros.h"

protected:
  int device_;
  int old_device_;
};
}
