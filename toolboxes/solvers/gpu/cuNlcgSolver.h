#pragma once

#include "nlcgSolver.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include "cuSolverUtils.h"

namespace Gadgetron{

  template <class T> class cuNlcgSolver : public nlcgSolver<cuNDArray<T> >
  {
  public:
    cuNlcgSolver() : nlcgSolver<cuNDArray<T> >() {}
    virtual ~cuNlcgSolver() {}
  };
}
