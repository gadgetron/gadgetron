#pragma once

#include "laplaceOperator.h"
#include "cuLinearOperator_macros.h"
#include "cuNDArray.h"
#include "solvers_export.h"

template <class REAL, class T, unsigned int D> 
class EXPORTSOLVERS cuLaplaceOperator : public laplaceOperator<REAL, D, cuNDArray<T> >
{
  
public:
  
  cuLaplaceOperator() : laplaceOperator< REAL, D, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuLaplaceOperator() {}

  virtual boost::shared_ptr< linearOperator< REAL, cuNDArray<T> > > clone(){
    return linearOperator< REAL, cuNDArray<T> >::clone(this);
  }

 protected:
  virtual int compute_laplace( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate );  

  DECLARE_LINEAR_OPERATOR_DEVICE_SUPPORT(cuLaplaceOperator)
};
