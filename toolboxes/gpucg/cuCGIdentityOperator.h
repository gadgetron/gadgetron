#pragma once

#include "cuCGMatrixOperator.h"
#include "ndarray_vector_td_utilities.h"

#include <cublas_v2.h>

template <class T> 
class cuCGIdentityOperator : public cuCGMatrixOperator<T>
{
 public:

 cuCGIdentityOperator( cublasHandle_t handle ) : handle_(handle) {}
  virtual ~cuCGIdentityOperator() {}
  
  virtual int mult_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false){

    if( accumulate ) 
      cuNDA_axpy( get_one<T>(), in, out, handle_ );
    else 
      *out = *in;

    return 0;
  }
  
  virtual int mult_MH(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false){

    if( accumulate ) 
      cuNDA_axpy( get_one<T>(), in, out, handle_ );
    else 
      *out = *in;

    return 0;
  }
  
  virtual int mult_MH_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false){

    if( accumulate ) 
      cuNDA_axpy( get_one<T>(), in, out, handle_ );
    else 
      *out = *in;
    
    return 0;
  }
  
 private:
  cublasHandle_t handle_;
};
