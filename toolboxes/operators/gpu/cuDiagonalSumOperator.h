/** \file cuDiagonalSumOperator.h
    \brief Sum of diagonal matrices, GPU instantiation.
*/

#pragma once

#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "diagonalSumOperator.h"

namespace Gadgetron{

  template <class T> class cuDiagonalSumOperator : public diagonalSumOperator< cuNDArray<T> >
  {
  public:
    cuDiagonalSumOperator() : diagonalSumOperator< cuNDArray<T> >() {}
    virtual ~cuDiagonalSumOperator() {}
  };
}
