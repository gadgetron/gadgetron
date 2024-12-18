/** \file cuPartialDerivativeOperator.h
    \brief Partial derivative regularization operator, GPU based.
*/

#pragma once

#include "cuNDArray_math.h"
#include "partialDerivativeOperator.h"

namespace Gadgetron{

  template <class T, unsigned int D> class cuPartialDerivativeOperator
    : public partialDerivativeOperator<D, cuNDArray<T> >
  {
  public:

    cuPartialDerivativeOperator() :
      partialDerivativeOperator< D, cuNDArray<T> >(0) {}

    cuPartialDerivativeOperator( size_t dimension ) :
      partialDerivativeOperator<D, cuNDArray<T> >( dimension ) {}

    virtual ~cuPartialDerivativeOperator() {}

    virtual void compute_partial_derivative( typename int64d<D>::Type stride, cuNDArray<T> *in,
                                             cuNDArray<T> *out, bool accumulate );

    virtual void compute_second_order_partial_derivative( typename int64d<D>::Type forwards_stride,
                                                          typename int64d<D>::Type adjoint_stride,
                                                          cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate );


  };
}
