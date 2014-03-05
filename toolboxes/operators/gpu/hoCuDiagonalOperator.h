/** \file hoCuDiagonalOperator.h
    \brief Diagonal matrix regularization operator for array type hoCuNDarray
*/

#pragma once

#include "hoCuNDArray_operators.h"
#include "hoCuNDArray_elemwise.h"
#include "hoCuNDArray_blas.h"
#include "diagonalOperator.h"

namespace Gadgetron{

  template <class T> class hoCuDiagonalOperator : public diagonalOperator< hoCuNDArray<T> >
  {
  public:
    hoCuDiagonalOperator() : diagonalOperator< hoCuNDArray<T> >() {}
    virtual ~hoCuDiagonalOperator() {}
  };
}
