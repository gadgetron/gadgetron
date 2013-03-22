#pragma once

#include "sbSolver.h"
#include "cuLwSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "encodingOperatorContainer.h"

template <class T> class cuSbLwSolver
  : public sbSolver<cuNDArray<typename realType<T>::type>, cuNDArray<T>, cuLwSolver<T> >
{
public:
  
  cuSbLwSolver() : sbSolver< cuNDArray<typename realType<T>::type>, cuNDArray<T>, cuLwSolver<T> >() {
    set_device(-1); 
  }

  virtual ~cuSbLwSolver() {}
  
#include "cuSbSolver_macros.h"
  
protected:
  int device_;
  int old_device_;
};
