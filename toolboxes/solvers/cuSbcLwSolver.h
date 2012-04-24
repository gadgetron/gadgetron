#pragma once

#include "sbcSolver.h"
#include "cuLwSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> class cuSbcLwSolver 
  : public sbcSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T>, cuLwSolver<REAL,T> >
{
public:
  
  cuSbcLwSolver() : sbcSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T>, cuLwSolver<REAL,T> >() { set_device(-1); }
  virtual ~cuSbcLwSolver() {}
  
#include "cuSbSolver_macros.h"

protected:
  int device_;
  int old_device_;
};
