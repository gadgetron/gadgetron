#pragma once

#include "cgPreconditioner.h"
#include "cuNDArray.h"
#include "gpusolvers_export.h"

#include <boost/shared_ptr.hpp>

namespace Gadgetron{

  template<class T> class EXPORTGPUSOLVERS cuCgPreconditioner : public cgPreconditioner< cuNDArray<T> >
  {
  public:
    
    cuCgPreconditioner() : cgPreconditioner< cuNDArray<T> >() {}
    virtual ~cuCgPreconditioner() {}
    
    virtual void set_weights( boost::shared_ptr< cuNDArray<T> > w );
    virtual void apply(cuNDArray<T>* in, cuNDArray<T>* out);
    
  protected:
    boost::shared_ptr< cuNDArray<T> > weights_;
  };
}
