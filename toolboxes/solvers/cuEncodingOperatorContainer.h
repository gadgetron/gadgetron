#pragma once

#include "encodingOperatorContainer.h"
#include "cuNDArray.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> class cuEncodingOperatorContainer 
  : public encodingOperatorContainer< REAL, cuNDArray<T> >
{
public:
  cuEncodingOperatorContainer() : encodingOperatorContainer< REAL, cuNDArray<T> >() {}
  virtual ~cuEncodingOperatorContainer() {}
  
  virtual bool operator_scal( REAL a, cuNDArray<T> *x ){
    return cuNDA_scal<REAL>(a,x);
  }
  
  virtual bool operator_axpy( REAL a,  cuNDArray<T> *x,  cuNDArray<T> *y ){
    return cuNDA_axpy<T>(a,x,y);
  }
  
  virtual boost::shared_ptr< linearOperator< REAL, cuNDArray<T> > > clone(){
    return linearOperator< REAL, cuNDArray<T> >::clone(this);
  }  
};
