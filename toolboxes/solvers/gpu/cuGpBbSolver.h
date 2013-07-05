#pragma once

#include "gpBbSolver.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"
#include "gpusolvers_export.h"

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/functional.h>

namespace Gadgetron{
  
  template <class T> class EXPORTGPUSOLVERS cuGpBbSolver : public gpBbSolver<cuNDArray<T> >
  {
  public:
    
    cuGpBbSolver() : gpBbSolver<cuNDArray<T> >() {}
    virtual ~cuGpBbSolver() {}
    
    virtual void solver_non_negativity_filter(cuNDArray<T> *x,cuNDArray<T> *g);    
  };
}
