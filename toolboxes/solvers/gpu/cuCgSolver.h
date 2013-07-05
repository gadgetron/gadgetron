/** \file cuCgSolver.h
    \brief Instantiation of the conjugate gradient solver on the cpu.

    The file cuCgSolver.h is a convienience wrapper for the device independent cgSolver class.
    The class cuCgSolver instantiates the cgSolver for the cuNDArray
    and the header otherwise includes other neccessary header files.
*/

#pragma once

#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "cgSolver.h"

namespace Gadgetron{
  
  /** \class cuCgSolver
      \brief Instantiation of the conjugate gradient solver on the cpu.
      
      The class cuCgSolver is a convienience wrapper for the device independent cgSolver class.
      cuCgSolver instantiates the cgSolver for type cuNDArray<T>.
  */
  template <class T> class cuCgSolver : public cgSolver< cuNDArray<T> >
  {
  public:    
    cuCgSolver() : cgSolver<cuNDArray<T> >() {}
    virtual ~cuCgSolver() {}
  };
}
