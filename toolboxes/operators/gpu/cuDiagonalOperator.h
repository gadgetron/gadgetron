/** \file cuDiagonalOperator.h
    \brief Diagonal matrix operator, GPU instantiation.
*/

#pragma once

#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "diagonalOperator.h"

namespace Gadgetron{

  template <class T> class cuDiagonalOperator : public diagonalOperator< cuNDArray<T> >
  {
  public:
    cuDiagonalOperator() : diagonalOperator< cuNDArray<T> >() {}
    virtual ~cuDiagonalOperator() {}
  };
}
