/** \file cuTv1DOperator.h
    \brief Total variation regularization operator, GPU based. Optimized 1D version.
*/

#pragma once

#include "cuNDArray_math.h"
#include "generalOperator.h"
#include "complext.h"

namespace Gadgetron{

  template<class T, unsigned int D> class cuTv1DOperator : public generalOperator< cuNDArray<T> >
  {


  public:

    typedef typename realType<T>::Type REAL;
    cuTv1DOperator() : generalOperator< cuNDArray<T> >(){
      limit_ = REAL(1e-8);
    }

    virtual ~cuTv1DOperator(){};

    void set_limit(REAL limit){
      limit_ = limit;
    }

    virtual void gradient(cuNDArray<T>*,cuNDArray<T>*, bool accumulate=false);
    virtual REAL magnitude(cuNDArray<T>*);

  protected:
    REAL limit_;
  };
}
