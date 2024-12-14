/** \file cuNonCartesianKtSenseOperator.h
    \brief Non-Cartesian kt-Sense operator, GPU based.
*/

#pragma once

#include "cuNonCartesianSenseOperator.h"

namespace Gadgetron{

  template<class REAL, unsigned int D>
  class cuNonCartesianKtSenseOperator : public cuNonCartesianSenseOperator<REAL,D>
  {

  public:

    typedef typename uint64d<D>::Type _uint64d;
    typedef typename reald<REAL,D>::Type _reald;

    cuNonCartesianKtSenseOperator() : cuNonCartesianSenseOperator<REAL,D>() {}
    virtual ~cuNonCartesianKtSenseOperator() {}

    virtual void mult_M( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate = false );
    virtual void mult_MH( cuNDArray< complext<REAL> >* in, cuNDArray< complext<REAL> >* out, bool accumulate = false );

  };
}
