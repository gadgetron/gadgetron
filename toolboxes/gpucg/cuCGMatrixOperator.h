#ifndef CUCGMATRIXOPERATOR_H
#define CUCGMATRIXOPERATOR_H

#include "cuNDArray.h"

template <class T> class cuCGMatrixOperator
{

 public:
  cuCGMatrixOperator() {}
  virtual ~cuCGMatrixOperator() {}

  virtual int mult_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
  virtual int mult_MH(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
  virtual int mult_MH_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
};

#endif //CUCGMATRIXOPERATOR_H
