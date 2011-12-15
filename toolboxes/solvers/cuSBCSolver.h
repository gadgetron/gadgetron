#pragma once

#include "sbcSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> class cuSBCSolver : public sbcSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >
{
public:
  
  cuSBCSolver( int device=-1 ) : sbcSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T> >() { set_device(device); }
  virtual ~cuSBCSolver() {}

#include "cuSBSolver_macros.h"

protected:
  int device_;
  int old_device_;
};
