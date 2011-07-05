#pragma once

#include "cuCGMatrixOperatorDevice.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> 
class cuCGIdentityOperator : public cuCGMatrixOperatorDevice<REAL,T>
{
 public:

  cuCGIdentityOperator( int device = -1 ) : cuCGMatrixOperatorDevice<REAL,T>(device) {}
  virtual ~cuCGIdentityOperator() {}
  
  virtual int mult_M_device(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false)
  {
    if( accumulate ) 
      cuNDA_axpy( get_one<T>(), in, out );
    else 
      *out = *in;

    return 0;
  }
  
  virtual int mult_MH_device(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false)
  {
    if( accumulate ) 
      cuNDA_axpy( get_one<T>(), in, out );
    else 
      *out = *in;
    
    return 0;
  }
  
  virtual int mult_MH_M_device(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false)
  {
    if( accumulate ) 
      cuNDA_axpy( get_one<T>(), in, out );
    else 
      *out = *in;
    
    return 0;
  }
};
