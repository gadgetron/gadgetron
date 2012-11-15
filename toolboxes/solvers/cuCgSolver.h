#pragma once

#include "cgSolver.h"
#include "cuNDArray.h"
#include "cuCgPreconditioner.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "cuGTBLAS.h"
#include <iostream>


template <class T> class cuCgSolver
	: public cgSolver<cuNDArray<T> >
{
public:

  cuCgSolver() : cgSolver<cuNDArray<T> >() {}
  virtual ~cuCgSolver() {}



    
protected:
  int device_;
  int old_device_;
  cuNDArray<T> *new_rhs;
};
