#ifndef CUCGPRECONDITIONER_H
#define CUCGPRECONDITIONER_H

#include "cuNDArray.h"
#include <boost/smart_ptr.hpp>

template <class T> class cuCGPreconditioner
{
 public:
  cuCGPreconditioner() {}
  virtual ~cuCGPreconditioner() {}
  virtual int apply(cuNDArray<T>* in, cuNDArray<T>* out) = 0;
};

template <class T> class cuCGPrecondWeight : public cuCGPreconditioner<T>
{
 public:
  virtual int set_weights( boost::shared_ptr< cuNDArray<T> > w ) {
    weights_ = w;
    return 0;
  }

  virtual int apply(cuNDArray<T>* in, cuNDArray<T>* out);

 protected:
  boost::shared_ptr< cuNDArray<T> > weights_;
};

#endif //CUCGPRECONDITIONER_H
