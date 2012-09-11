#pragma once

#include "diagonalOperator.h"
#include "cuLinearOperator_macros.h"
#include "cuNDArray.h"

template <class REAL, class T> class cuDiagonalOperator 
  : public diagonalOperator< REAL, cuNDArray<T> >
{

public:

  cuDiagonalOperator() : diagonalOperator< REAL, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuDiagonalOperator() {}
     
  virtual int mult_MH_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false )
  {
    _set_device();
    int res = diagonalOperator< REAL, cuNDArray<T> >::mult_MH_M( in, out, accumulate );
    _restore_device();
    return res;
  }
    
  virtual bool operator_clear( cuNDArray<T> *x )
  {
    return cuNDA_clear( x );
  }
 
  virtual bool operator_axpy( cuNDArray<T> *a,  cuNDArray<T> *x,  cuNDArray<T> *y )
  {
    return cuNDA_axpy( a, x, y );
  }

  virtual boost::shared_ptr< linearOperator< T, cuNDArray<T> > > clone()
  {
    return linearOperator< REAL, cuNDArray<T> >::clone(this);
  }
  
  DECLARE_LINEAR_OPERATOR_DEVICE_SUPPORT(cuDiagonalOperator)
};
