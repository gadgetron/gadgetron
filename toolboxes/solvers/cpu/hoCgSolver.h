/** \file hoCgSolver.h
    \brief Instantiation of the conjugate gradient solver on the cpu.

    The file hoCgSolver.h is a convienience wrapper for the device independent cgSolver class.
    The class hoCgSolver instantiates the cgSolver for the hoNDArray
    and the header otherwise includes other neccessary header files.
*/

#pragma once

#include "cgSolver.h"
#include "hoNDArray_math.h"

namespace Gadgetron{

  /** \class hoCgSolver
      \brief Instantiation of the conjugate gradient solver on the cpu.
      
      The class hoCgSolver is a convienience wrapper for the device independent cgSolver class.
      hoCgSolver instantiates the cgSolver for type hoNDArray<T>.
  */
  template <class T> class hoCgSolver : public cgSolver< hoNDArray<T> >
  {
  public:
    hoCgSolver() : cgSolver<hoNDArray<T> >() {}
    virtual ~hoCgSolver() {}
  };
}
