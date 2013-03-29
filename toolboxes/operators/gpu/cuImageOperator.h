#pragma once

#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "complext.h"
#include "imageOperator.h"

namespace Gadgetron{

  template <class T> class cuImageOperator : public imageOperator< cuNDArray<typename realType<T>::Type >, cuNDArray<T> >
  {
  public:
    cuImageOperator() : imageOperator< cuNDArray<typename realType<T>::Type >, cuNDArray<T> >() {}
    virtual ~cuImageOperator() {}
  };
}
