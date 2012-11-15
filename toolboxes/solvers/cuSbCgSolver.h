#pragma once

#include "sbSolver.h"
#include "cuCgSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "complext.h"

template <class T> class cuSbCgSolver
  : public sbSolver<cuNDArray<typename realType<T>::type >, cuNDArray<T>, cuCgSolver<T> >
{
public:
  
  cuSbCgSolver() : sbSolver<cuNDArray<typename realType<T>::type >, cuNDArray<T>, cuCgSolver<T> >() {
  }
  
  virtual ~cuSbCgSolver() {}
};
