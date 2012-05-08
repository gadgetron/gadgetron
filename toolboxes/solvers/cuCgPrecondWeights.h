#pragma once

#include "cuCgPreconditioner.h"
#include <boost/smart_ptr.hpp>

template<class T> class EXPORTSOLVERS cuCgPrecondWeights : public cuCgPreconditioner<T>
{
 public:

  cuCgPrecondWeights( int device = -1 ) : cuCgPreconditioner<T>(device) {}
  virtual ~cuCgPrecondWeights() {}
  
  virtual int set_weights( boost::shared_ptr< cuNDArray<T> > w );
  virtual int apply(cuNDArray<T>* in, cuNDArray<T>* out);

 protected:
  boost::shared_ptr< cuNDArray<T> > weights_;
};
