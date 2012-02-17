#pragma once

#include "laplaceOperator.h"
#include "cuMatrixOperator_macros.h"
#include "cuNDArray.h"

template <class REAL, class T, unsigned int D> 
class EXPORTSOLVERS cuLaplaceOperator : public laplaceOperator<REAL, D, cuNDArray<T> >
{
  
public:
  
  cuLaplaceOperator() : laplaceOperator< REAL, D, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuLaplaceOperator() {}

 protected:
  virtual int compute_laplace( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate );  

  DECLARE_MATRIX_OPERATOR_DEVICE_SUPPORT(cuLaplaceOperator)
};
