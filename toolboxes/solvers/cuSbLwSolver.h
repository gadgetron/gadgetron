#pragma once

#include "sbSolver.h"
#include "cuLwSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "cuEncodingOperatorContainer.h"

template <class REAL, class T> class cuSbLwSolver 
  : public sbSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T>, cuLwSolver<REAL,T>, cuEncodingOperatorContainer<REAL,T> >
{
public:
  
  cuSbLwSolver() : sbSolver< REAL, T, cuNDArray<REAL>, cuNDArray<T>, cuLwSolver<REAL,T>, cuEncodingOperatorContainer<REAL,T> >() { 
    set_device(-1); 
  }

  virtual ~cuSbLwSolver() {}
  
#include "cuSbSolver_macros.h"
  
protected:
  int device_;
  int old_device_;
};
