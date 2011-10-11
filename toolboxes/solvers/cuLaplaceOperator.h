#pragma once

#include "cuNDArray.h"
#include "laplaceOperator.h"

template <class REAL, class T, unsigned int D> class EXPORTSOLVERS cuLaplaceOperator : public laplaceOperator<REAL, D, cuNDArray<T> >
{
  
public:
  
  cuLaplaceOperator() : laplaceOperator< REAL, D, cuNDArray<T> >() {}
  virtual ~cuLaplaceOperator() {}
 protected:
  virtual int compute_laplace( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate );  


};
