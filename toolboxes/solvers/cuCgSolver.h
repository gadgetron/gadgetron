#pragma once

#include "cuNDArray.h"
#include "cgSolver.h"
#include "cuGTBLAS.h"
namespace Gadgetron{
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
}
