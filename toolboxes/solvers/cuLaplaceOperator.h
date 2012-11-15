#pragma once

#include "laplaceOperator.h"
#include "cuLinearOperator_macros.h"
#include "cuNDArray.h"
#include "solvers_export.h"

template < class T, unsigned int D>
class EXPORTSOLVERS cuLaplaceOperator : public laplaceOperator<D, cuNDArray<T> >
{
  
public:
  
  cuLaplaceOperator() : laplaceOperator< D, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuLaplaceOperator() {}

  virtual boost::shared_ptr< linearOperator< cuNDArray<T> > > clone(){
    return linearOperator<cuNDArray<T> >::clone(this);
  }

 protected:
  virtual void compute_laplace( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate );

  DECLARE_LINEAR_OPERATOR_DEVICE_SUPPORT(cuLaplaceOperator)
};
