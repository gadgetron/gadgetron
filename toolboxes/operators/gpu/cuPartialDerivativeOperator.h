#pragma once

#include "partialDerivativeOperator.h"
#include "cuNDArray.h"
#include "gpuoperators_export.h"

namespace Gadgetron{

  template <class T, unsigned int D> class EXPORTGPUOPERATORS cuPartialDerivativeOperator 
    : public partialDerivativeOperator<D, cuNDArray<T> >
  {
  public:
    
    cuPartialDerivativeOperator() : 
      partialDerivativeOperator< D, cuNDArray<T> >(0) {}
    
    cuPartialDerivativeOperator( unsigned int dimension ) : 
      partialDerivativeOperator<D, cuNDArray<T> >( dimension ) {}
    
    virtual ~cuPartialDerivativeOperator() {}
    
    virtual void compute_partial_derivative( typename intd<D>::Type stride, cuNDArray<T> *in,
					     cuNDArray<T> *out, bool accumulate );  

    virtual void compute_second_order_partial_derivative( typename intd<D>::Type forwards_stride,
							  typename intd<D>::Type adjoint_stride, 
							  cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate );  

    virtual boost::shared_ptr< linearOperator< cuNDArray<T> > > clone() {
      return linearOperator< cuNDArray<T> >::clone(this);
    }    
  };
}
