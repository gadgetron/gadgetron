#pragma once

#include "cgSolver.h"

#include "cgSolver.h"
#include "hoNDArray_operators.h"
#include "hoNDArray_elemwise.h"
#include "hoCuNDArray_blas.h"

namespace Gadgetron{

  /** \class hoCuCgSolver
      \brief Instantiation of the conjugate gradient solver on the cpu.

      The class hoCuCgSolver is a convienience wrapper for the device independent cgSolver class.
      hoCuCgSolver instantiates the cgSolver for type hoNDArray<T>.
  */
  template <class T> class hoCuCgSolver : public cgSolver< hoCuNDArray<T> >
  {
  public:
    hoCuCgSolver() : cgSolver<hoCuNDArray<T> >(), _it(0) {}
    virtual ~hoCuCgSolver() {}

    virtual void solver_dump(hoCuNDArray<T>* x){
    	std::stringstream ss;
			ss << "iteration-" << _it << ".real";
			write_nd_array(x,ss.str().c_str());
			_it++;
    }

  private:
    int _it;
  };
}
