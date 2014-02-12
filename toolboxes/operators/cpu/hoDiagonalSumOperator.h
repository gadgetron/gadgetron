/** \file hoDiagonalSumOperator.h
    \brief Sum of diagonal matrices operator, CPU instantiation.
*/

#pragma once

#include "hoNDArray_operators.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_blas.h"
#include "diagonalSumOperator.h"

namespace Gadgetron{

  template <class T> class hoDiagonalSumOperator : public diagonalSumOperator< hoNDArray<T> >
  {
  public:
    hoDiagonalSumOperator() : diagonalSumOperator< hoNDArray<T> >() {}
    virtual ~hoDiagonalSumOperator() {}
  };
}
