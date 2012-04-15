#pragma once

#include "sbSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> class cuSbSolver 
	: public sbSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >
{
public:
  
  cuSbSolver() : sbSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuSbSolver() {}

#include "cuSbSolver_macros.h"

protected:
  int device_;
  int old_device_;
};
