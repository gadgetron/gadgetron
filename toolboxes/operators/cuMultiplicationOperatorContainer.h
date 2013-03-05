#pragma once

#include "multiplicationOperatorContainer.h"
#include "cuNDArray.h"

namespace Gadgetron{
template <class REAL, class T> class cuMultiplicationOperatorContainer 
  : public multiplicationOperatorContainer< REAL, cuNDArray<T> >
{
public:
  cuMultiplicationOperatorContainer() : multiplicationOperatorContainer< REAL, cuNDArray<T> >() {}
  virtual ~cuMultiplicationOperatorContainer() {}
  
  virtual boost::shared_ptr< linearOperator< REAL, cuNDArray<T> > > clone(){
    return linearOperator< REAL, cuNDArray<T> >::clone(this);
  }  
};
}
