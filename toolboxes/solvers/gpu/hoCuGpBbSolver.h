#pragma once

#include "gpBbSolver.h"
#include "hoNDArray_math.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "hoSolverUtils.h"

#ifdef USE_OMP
#include <omp.h>
#endif

namespace Gadgetron{

  template <class T> class hoCuGpBbSolver : public gpBbSolver< hoCuNDArray<T> >
  {  
  public:

    hoCuGpBbSolver() : gpBbSolver< hoCuNDArray<T> >() {};
    virtual ~hoCuGpBbSolver() {};

  };
}
