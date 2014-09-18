#pragma once

#include "gpBbSolver.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"


#include "cuSolverUtils.h"

namespace Gadgetron{

  template <class T> class cuGpBbSolver : public gpBbSolver<cuNDArray<T> >
  {
  public:

    cuGpBbSolver() : gpBbSolver<cuNDArray<T> >() {}
    virtual ~cuGpBbSolver() {}
  };
}
