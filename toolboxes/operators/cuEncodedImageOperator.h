#pragma once

#include "encodedImageOperator.h"
#include "cuLinearOperator_macros.h"
#include "cuNDArray.h"

namespace Gadgetron{
template <class T> class cuEncodedImageOperator
	: public encodedImageOperator<cuNDArray<typename realType<T>::type >, cuNDArray<T> >
{
		typedef typename realType<T>::type REAL;
 public:

  cuEncodedImageOperator() : encodedImageOperator<cuNDArray<REAL>, cuNDArray<T> >() { set_device(-1); }
  virtual ~cuEncodedImageOperator() {}

  virtual void compute( cuNDArray<T> *image )
  {
    _set_device();
    encodedImageOperator<cuNDArray<REAL>, cuNDArray<T> >::compute( image );
    _restore_device();
  }

  virtual void mult_MH_M( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false )
  {
    _set_device();
    encodedImageOperator<cuNDArray<REAL>, cuNDArray<T> >::mult_MH_M( in, out, accumulate );
    _restore_device();
  }
  
  
  virtual boost::shared_ptr< linearOperator<cuNDArray<T> > > clone()
  {
    return linearOperator<cuNDArray<T> >::clone(this);
  }
  
  DECLARE_LINEAR_OPERATOR_DEVICE_SUPPORT(cuEncodedImageOperator)
};
}
