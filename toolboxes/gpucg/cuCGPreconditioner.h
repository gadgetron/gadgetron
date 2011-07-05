#ifndef CUCGPRECONDITIONER_H
#define CUCGPRECONDITIONER_H

#pragma once
#include "gadgetron_export.h"
#include "cuNDArray.h"

#include <boost/smart_ptr.hpp>

template <class T> class EXPORTGPUCG cuCGPreconditioner
{
 public:
  cuCGPreconditioner() {}
  virtual ~cuCGPreconditioner() {}
  virtual int apply(cuNDArray<T>* in, cuNDArray<T>* out) = 0;

  void* operator new (size_t bytes) { return ::new char[bytes]; }
  void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
  void * operator new(size_t s, void * p) { return p; }
};

template <class T> class EXPORTGPUCG cuCGPrecondWeight : public cuCGPreconditioner<T>
{
 public:
  cuCGPrecondWeight( int device = -1 )
  {
    if( device<0 ){
      if( cudaGetDevice( &device_ ) != cudaSuccess ){
	std::cerr << "cuCGMatrixOperatorDevice: unable to get current device." << std::endl ;
	device_ = 0;
      }      
    }
    else
      device_ = device;    
  }

  virtual ~cuCGPrecondWeight() {}

  virtual int set_weights( boost::shared_ptr< cuNDArray<T> > w ) {
    weights_ = w;
    return 0;
  }

  virtual int apply(cuNDArray<T>* in, cuNDArray<T>* out);

 protected:
  boost::shared_ptr< cuNDArray<T> > weights_;
  int device_;
};

#endif //CUCGPRECONDITIONER_H
