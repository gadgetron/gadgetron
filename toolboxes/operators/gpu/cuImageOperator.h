#pragma once

#include "imageOperator.h"
#include "cuNDArray.h"

namespace Gadgetron{
  template <class T> class cuImageOperator : public imageOperator<cuNDArray<typename realType<T>::Type>, cuNDArray<T> >
  {    
  public:    
    cuImageOperator() : imageOperator< cuNDArray<typename realType<T>::Type >, cuNDArray<T> >() {}
    virtual ~cuImageOperator() {}     
  };
}
