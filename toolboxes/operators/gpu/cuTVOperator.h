#pragma once

#include "generalOperator.h"
#include "cuNDArray.h"
#include "complext.h"
#include "gpuoperators_export.h"

namespace Gadgetron{

  template<class T, unsigned int D> class EXPORTGPUOPERATORS cuTVOperator : public generalOperator<cuNDArray<T> > 
  {
  protected:
    typedef typename realType<T>::Type REAL;
    
  public:

    cuTVOperator() : generalOperator<cuNDArray<T> >(){
      limit_ = REAL(1e-8);
    }

    virtual ~cuTVOperator(){};

    void set_limit(REAL limit){
      limit_ = limit;
    }

    virtual void gradient(cuNDArray<T>*,cuNDArray<T>*, bool accumulate=false);
  protected:

  protected:    
    REAL limit_;
  };
}
