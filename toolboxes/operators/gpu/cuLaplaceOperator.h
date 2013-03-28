#pragma once

#include "laplaceOperator.h"
#include "cuNDArray.h"
#include "gpuoperators_export.h"

namespace Gadgetron{

  template < class T, unsigned int D> class EXPORTGPUOPERATORS cuLaplaceOperator : public laplaceOperator<D, cuNDArray<T> >
  {    
  public:
    
    cuLaplaceOperator() : laplaceOperator< D, cuNDArray<T> >() {}
    virtual ~cuLaplaceOperator() {}
    
    virtual boost::shared_ptr< linearOperator< cuNDArray<T> > > clone(){
      return linearOperator<cuNDArray<T> >::clone(this);
    }
    
  protected:
    virtual void compute_laplace( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate );    
  };
}
