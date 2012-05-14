#pragma once

#include "imageOperator.h"
#include "cuLinearOperator_macros.h"
#include "cuNDArray.h"

template <class REAL, class T> class cuImageOperator 
	: public imageOperator< REAL, cuNDArray<REAL>, cuNDArray<T> >
{

 public:

  cuImageOperator() : imageOperator< REAL, cuNDArray<REAL>, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuImageOperator() {}
     
  virtual int compute( cuNDArray<T> *image )
  {
    _set_device();
    int res = imageOperator< REAL, cuNDArray<REAL>, cuNDArray<T> >::compute( image );
    _restore_device();
    return res;
  }

  virtual int mult_MH_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false )
  {
    _set_device();
    int res = imageOperator< REAL, cuNDArray<REAL>, cuNDArray<T> >::mult_MH_M( in, out, accumulate );
    _restore_device();
    return res;
  }
  
  virtual void operator_scal( REAL scale, cuNDArray<T> *x )
  {
    cuNDA_scal<T>( scale*T(1), x );
  }

  virtual void operator_reciprocal( cuNDArray<REAL> *x )
  {
    cuNDA_reciprocal<REAL>(x);
  }

  virtual REAL operator_asum( cuNDArray<T> *x )
  {
    return cuNDA_asum<REAL,T>(x);
  }

  virtual boost::shared_ptr< cuNDArray<REAL> > operator_abs( cuNDArray<T> *x )
  {
	  boost::shared_ptr< cuNDArray<REAL> > res =cuNDA_cNorm<REAL,T>( x, CUNDA_NDARRAY_DEVICE, CUNDA_NDARRAY_DEVICE);
	  cuNDA_sqrt(res.get());
	  return res;
  }
  
  virtual bool operator_clear(  cuNDArray<T> *x )
  {
    return cuNDA_clear(x);
  }
 
  virtual bool operator_axpy( cuNDArray<REAL> *a,  cuNDArray<T> *x,  cuNDArray<T> *y )
  {
    return cuNDA_axpy( a, x, y );
  }

  virtual boost::shared_ptr< linearOperator< REAL, cuNDArray<T> > > clone()
  {
    return linearOperator< REAL, cuNDArray<T> >::clone(this);
  }
  
  DECLARE_LINEAR_OPERATOR_DEVICE_SUPPORT(cuImageOperator)
};
