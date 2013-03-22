#pragma once

#include "imageOperator.h"
#include "cuLinearOperator_macros.h"
#include "cuNDArray.h"
#include "ndarray_vector_td_utilities.h"

namespace Gadgetron{
template <class T> class cuImageOperator
  : public imageOperator<cuNDArray<typename realType<T>::type>, cuNDArray<T> >
{

public:

  cuImageOperator() : imageOperator< cuNDArray<typename realType<T>::type >, cuNDArray<T> >() {}
  virtual ~cuImageOperator() {}
     
  
};
}
