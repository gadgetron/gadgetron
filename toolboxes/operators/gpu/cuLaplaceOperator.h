/** \file cuLaplaceOperator.h
    \brief Laplace regularization operator, GPU based.
*/

#pragma once

#include "gpuoperators_export.h"
#include "cuNDArray_math.h"
#include "laplaceOperator.h"

namespace Gadgetron{

  template < class T, unsigned int D> class EXPORTGPUOPERATORS cuLaplaceOperator : public laplaceOperator<D, cuNDArray<T> >
  {    
  public:
    
    cuLaplaceOperator() : laplaceOperator< D, cuNDArray<T> >() {}
    virtual ~cuLaplaceOperator() {}
    

  protected:
    virtual void compute_laplace( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate );    
  };
}
