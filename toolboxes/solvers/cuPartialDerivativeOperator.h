#pragma once

#include "partialDerivativeOperator.h"
#include "cuLinearOperator_macros.h"
#include "cuNDArray.h"

template <class T, unsigned int D> class EXPORTSOLVERS cuPartialDerivativeOperator
	: public partialDerivativeOperator<D, cuNDArray<T> >
{
 public:
  
  cuPartialDerivativeOperator() : 
    partialDerivativeOperator< D, cuNDArray<T> >(0) { set_device(-1); }

  cuPartialDerivativeOperator( unsigned int dimension ) : 
    partialDerivativeOperator<D, cuNDArray<T> >( dimension ) { set_device(-1); }
  
  virtual ~cuPartialDerivativeOperator() {}

  virtual boost::shared_ptr< linearOperator< cuNDArray<T> > > clone() {
    return linearOperator<cuNDArray<T> >::clone(this);
  }
  
  virtual void compute_partial_derivative( typename intd<D>::Type stride, cuNDArray<T> *in,
					  cuNDArray<T> *out, bool accumulate );  

  virtual void compute_second_order_partial_derivative( typename intd<D>::Type forwards_stride,
						       typename intd<D>::Type adjoint_stride, 
						       cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate );  
  
  DECLARE_LINEAR_OPERATOR_DEVICE_SUPPORT(cuPartialDerivativeOperator)
};
