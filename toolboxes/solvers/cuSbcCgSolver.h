#pragma once

#include "sbcSolver.h"
#include "cuCgSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"
#include "complext.h"

template <class T> class cuSbcCgSolver
  : public sbcSolver<
		      cuNDArray<typename realType<T>::type >,
		      cuNDArray<T>,
		      cuCgSolver<T> >
{
public:
  
  cuSbcCgSolver() : sbcSolver<cuNDArray<typename realType<T>::type >, cuNDArray<T>,
			       cuCgSolver<T> >() {}
  virtual ~cuSbcCgSolver() {}

};
