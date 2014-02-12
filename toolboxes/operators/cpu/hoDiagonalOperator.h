/** \file hoDiagonalOperator.h
    \brief Diagonal matrix operator, CPU instantiation.
*/

#pragma once

#include "hoNDArray_operators.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_blas.h"
#include "diagonalOperator.h"

namespace Gadgetron{

  template <class T> class hoDiagonalOperator : public diagonalOperator< hoNDArray<T> >
  {
  public:
    hoDiagonalOperator() : diagonalOperator< hoNDArray<T> >() {}
    virtual ~hoDiagonalOperator() {}
  };
}
