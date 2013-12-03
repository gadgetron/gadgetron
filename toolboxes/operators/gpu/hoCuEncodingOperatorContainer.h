#pragma once

#include "hoCuNDArray.h"
#include "hoCuNDArray_operators.h"
#include "hoCuNDArray_elemwise.h"
#include "hoCuNDArray_blas.h"
#include "encodingOperatorContainer.h"

namespace Gadgetron{
  
  template<class T> class hoCuEncodingOperatorContainer 
    : public encodingOperatorContainer< hoCuNDArray<T> >
  {
  public:
    hoCuEncodingOperatorContainer() : encodingOperatorContainer< hoCuNDArray<T> >() {}
    virtual ~hoCuEncodingOperatorContainer() {}
    
    virtual boost::shared_ptr< linearOperator< hoCuNDArray<T> > > clone(){
      return linearOperator< hoCuNDArray<T> >::clone(this);
    }  
  }; 
}
