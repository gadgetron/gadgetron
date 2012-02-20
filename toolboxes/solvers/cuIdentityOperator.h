#pragma once

#include "identityOperator.h"
#include "cuMatrixOperator_macros.h"
#include "ndarray_vector_td_utilities.h"

template <class REAL, class T> 
class cuIdentityOperator : public identityOperator< REAL, cuNDArray<T> >
{
 public:

  cuIdentityOperator() : identityOperator< REAL, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuIdentityOperator() {}
  
  virtual bool operator_xpy( cuNDArray<T> *x, cuNDArray<T> *y )
  { 
    int ret1 = _set_device();
    bool ret2;
    if( ret1 == 0 )
      ret2 = cuNDA_axpy( T(1), x, y, CUNDA_CURRENT_DEVICE );
    else 
      ret2 = false;
    ret1 = _restore_device();

    if( ret1 == 0 && ret2 )      
      return true;
    else{
      std::cout << std::endl << "Error :: cuIdentityOperator :: operator_xpy failed" << std::endl;
      return false;
    }
  }

  DECLARE_MATRIX_OPERATOR_DEVICE_SUPPORT(cuIdentityOperator);
};
