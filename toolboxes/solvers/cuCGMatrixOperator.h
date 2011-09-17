#ifndef CUCGMATRIXOPERATOR_H
#define CUCGMATRIXOPERATOR_H
#pragma once

#include "cgMatrixOperator.h"
#include "cuNDArray.h"

template <class REAL, class T> class cuCGMatrixOperator : public cgMatrixOperator< REAL, cuNDArray<T> >
{
 public:

  cuCGMatrixOperator( int device = -1 ) : cgMatrixOperator< REAL, cuNDArray<T> >() { set_device(device); }  
  virtual ~cuCGMatrixOperator() {}

  virtual void set_device( int device )
  {
    if( device<0 ){
      if( cudaGetDevice( &device_ ) != cudaSuccess ){
	std::cerr << "cuCGMatrixOperator::set_device: unable to get current device." << std::endl;
	device_ = 0;
      }      
    }
    else
      device_ = device;
  }

protected:
  virtual int set_device()
  {
    if( cudaGetDevice( &old_device_ ) != cudaSuccess ){
      std::cerr << "cuCGMatrixOperator::set_device: unable to get current device." << std::endl ;
      return -1;
    }
    if( device_ != old_device_ && cudaSetDevice(device_) != cudaSuccess) {
      std::cerr << "cuCGMatrixOperator::set_device: unable to set device " << device_ << std::endl;
      return -1;
    }
    return 0;
  }

  virtual int restore_device()
  {
    if( device_ != old_device_ && cudaSetDevice(old_device_) != cudaSuccess) {
      std::cerr << "cuCGMatrixOperator::restore_device: unable to set device " << old_device_ << std::endl;
      return -1;
    }
    return 0;
  }

protected:
  int device_;

private:
  int old_device_;

};

#endif //CUCGMATRIXOPERATOR_H
