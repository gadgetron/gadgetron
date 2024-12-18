/** \file cuSenseOperator.h
    \brief Base class for the GPU based Sense operators
*/

#pragma once

#include "senseOperator.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "vector_td.h"
#include "complext.h"

namespace Gadgetron{

  template<class REAL, unsigned int D> class cuSenseOperator : public senseOperator< cuNDArray< complext<REAL> >, D >
  {

  public:

    cuSenseOperator() : senseOperator<cuNDArray< complext<REAL> >,D >() {}
    virtual ~cuSenseOperator() {}

    virtual void mult_M( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false ) = 0;
    virtual void mult_MH( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out, bool accumulate = false ) = 0;

    virtual void mult_csm( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out );
    virtual void mult_csm_conj_sum( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out );
  };
}
