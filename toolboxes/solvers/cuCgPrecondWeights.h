#pragma once

#include "cuCgPreconditioner.h"
#include <boost/smart_ptr.hpp>
namespace Gadgetron{
template<class T> class EXPORTSOLVERS cuCgPrecondWeights : public cuCgPreconditioner<T>
{
 public:

  cuCgPrecondWeights( int device = -1 ) : cuCgPreconditioner<T>(device) {}
  virtual ~cuCgPrecondWeights() {}
  
  virtual void set_weights( boost::shared_ptr< cuNDArray<T> > w );
  virtual void apply(cuNDArray<T>* in, cuNDArray<T>* out);

 protected:
  boost::shared_ptr< cuNDArray<T> > weights_;
};
}
