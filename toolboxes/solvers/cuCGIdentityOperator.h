#pragma once

#include "cuCGMatrixOperator.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> 
class cuCGIdentityOperator : public cuCGMatrixOperator<REAL,T>
{
 public:

  cuCGIdentityOperator( int device = -1 ) : cuCGMatrixOperator<REAL,T>(device) {}
  virtual ~cuCGIdentityOperator() {}
  
  virtual int mult_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false)
  {
    int ret1 = 0;
    bool ret2 = true;    
    
    if( accumulate ){
      ret1 = this->set_device();
      if( ret1 == 0 )
	ret2 = cuNDA_axpy( get_one<T>(), in, out, CUNDA_CURRENT_DEVICE );
      else 
	ret2 = false;
      ret1 = this->restore_device();
    }
    else 
      *out = *in;
    
    if( ret1 == 0 && ret2 )      
      return 0;
    else{
      std::cout << std::endl << "cuCGIdentityOperator::mult failed" << std::endl;
      return -1;
    }
  }
  
  virtual int mult_MH(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false)
  {
    return mult_M(in, out, accumulate);
  }
  
  virtual int mult_MH_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false)
  {
    return mult_M(in, out, accumulate);
  }
};
