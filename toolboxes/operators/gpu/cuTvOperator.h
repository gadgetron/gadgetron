/** \file cuTvOperator.h
    \brief Total variation regularization operator, GPU based.
*/

#pragma once

#include "cuNDArray_math.h"
#include "generalOperator.h"

#include "complext.h"

namespace Gadgetron{

  template<class T, unsigned int D> class cuTvOperator
    : public generalOperator<cuNDArray<T> >
  {


  public:
    typedef typename realType<T>::Type REAL;

    cuTvOperator() : generalOperator<cuNDArray<T> >(){
      limit_ = REAL(0);
    }

    virtual ~cuTvOperator(){};

    void set_limit(REAL limit){
      limit_ = limit;
    }

    virtual void gradient(cuNDArray<T>*,cuNDArray<T>*, bool accumulate=false);
    virtual REAL magnitude(cuNDArray<T>*);

  protected:

  protected:
    REAL limit_;
  };
}
