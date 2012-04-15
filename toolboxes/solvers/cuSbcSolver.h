#pragma once

#include "sbcSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> class cuSbcSolver 
	: public sbcSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >
{
public:
  
  cuSbcSolver() : sbcSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuSbcSolver() {}

#include "cuSbSolver_macros.h"

protected:
  int device_;
  int old_device_;
};
