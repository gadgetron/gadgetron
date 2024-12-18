#pragma once

#include "linearOperator.h"
#include "vector_td.h"
#include "cuNDArray.h"


namespace Gadgetron {
template <class T, unsigned int D> class cuPartialDerivativeOperator2
	: public linearOperator<cuNDArray<T> >
{

public:

  cuPartialDerivativeOperator2() : linearOperator<cuNDArray<T> >() {}
  virtual ~cuPartialDerivativeOperator2() {}

  virtual void mult_M( cuNDArray<T> *in,cuNDArray<T> *out, bool accumulate = false );


  virtual void mult_MH(cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate = false );

  virtual void mult_MH_M( cuNDArray<T> *in, cuNDArray<T>*out, bool accumulate = false );

};
}
