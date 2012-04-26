#pragma once

#include "encodingOperatorContainer.h"
#include "cuNDArray.h"

template <class REAL, class T> class cuEncodingOperatorContainer 
  : public encodingOperatorContainer< REAL, cuNDArray<T> >
{
public:
  cuEncodingOperatorContainer() : encodingOperatorContainer< REAL, cuNDArray<T> >() {}
  virtual ~cuEncodingOperatorContainer() {}
  
  virtual boost::shared_ptr< linearOperator< REAL, cuNDArray<T> > > clone()
  {
    return linearOperator< REAL, cuNDArray<T> >::clone(this);
  }  
};
