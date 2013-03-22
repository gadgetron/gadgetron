#ifndef CUCGPRECONDITIONER_H
#define CUCGPRECONDITIONER_H
#pragma once

#include "cgPreconditioner.h"
#include "cuNDArray.h"
namespace Gadgetron{


template <class T> class cuCgPreconditioner 
	: public cgPreconditioner< cuNDArray<T> >
{
 public:
  cuCgPreconditioner( int device = -1 ) : cgPreconditioner< cuNDArray<T> >()
  {
    if( device<0 ){
      if( cudaGetDevice( &device_ ) != cudaSuccess ){
	std::cerr << "cuCgPreconditioner: unable to get current device." << std::endl ;
	device_ = 0;
      }
    }
    else
      device_ = device;
  }

  virtual ~cuCgPreconditioner() {}

protected:
  virtual int set_device()
  {
    if( cudaGetDevice( &old_device_ ) != cudaSuccess ){
      std::cerr << "cuCgPreconditioner::set_device: unable to get current device." << std::endl ;
      return -1;
    }
    if( device_ != old_device_ && cudaSetDevice(device_) != cudaSuccess) {
      std::cerr << "cuCgPreconditioner::set_device: unable to set device " << device_ << std::endl;
      return -1;
    }
    return 0;
  }

  virtual int restore_device()
  {
    if( device_ != old_device_ && cudaSetDevice(old_device_) != cudaSuccess) {
      std::cerr << "cuCgPreconditioner::restore_device: unable to set device " << old_device_ << std::endl;
      return -1;
    }
    return 0;
  }

 protected:
  int device_;
  
 private:
  int old_device_;

};

}
#endif //CUCGPRECONDITIONER_H
