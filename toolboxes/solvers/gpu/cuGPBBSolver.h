#pragma once

#include "ndarray_vector_td_utilities.h"
#include "cuGTBLAS.h"
#include "gpBBSolver.h"
#include "cuNDArray.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/functional.h>

namespace Gadgetron{
template <class T> class cuGPBBSolver
: public gpBBSolver<cuNDArray<T> >{

public:
	typedef typename realType<T>::type REAL;
	cuGPBBSolver() : gpBBSolver<cuNDArray<T> >() { };
	virtual	~cuGPBBSolver(){};


	  virtual void solver_non_negativity_filter(cuNDArray<T> *x,cuNDArray<T> *g);

	  virtual void solver_reciprocal_clamp( cuNDArray<T>* x,REAL threshold) ;




};
}
