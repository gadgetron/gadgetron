#pragma once

#include "sbSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> class cuSBSolver 
	: public sbSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >
{
public:
  
  cuSBSolver( int device=-1 ) : sbSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >() { set_device(device); }
  virtual ~cuSBSolver() {}

#include "cuSBSolver_macros.h"

protected:
  int device_;
  int old_device_;
};
