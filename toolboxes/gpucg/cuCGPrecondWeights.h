#pragma once

#include "cuCGPreconditioner.h"
#include <boost/smart_ptr.hpp>
#include "gadgetron_export.h"

template<class T> class EXPORTGPUCG cuCGPrecondWeights : public cuCGPreconditioner<T>
{
 public:

  cuCGPrecondWeights( int device = -1 ) : cuCGPreconditioner<T>(device) {}
  virtual ~cuCGPrecondWeights() {}
  
  virtual int set_weights( boost::shared_ptr< cuNDArray<T> > w );
  virtual int apply(cuNDArray<T>* in, cuNDArray<T>* out);

 protected:
  boost::shared_ptr< cuNDArray<T> > weights_;
};
