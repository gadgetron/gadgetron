#ifndef CUCGMATRIXOPERATOR_H
#define CUCGMATRIXOPERATOR_H
#pragma once

#include "cuNDArray.h"
#include "vector_td_utilities.h"

template <class REAL, class T> class cuCGMatrixOperator
{
 public:

  cuCGMatrixOperator( int device = -1 ) 
  { 
    weight_ = get_one<REAL>(); 
    set_device(device);
  }
  
  virtual ~cuCGMatrixOperator() {}

  inline void set_weight( REAL weight ){ weight_ = weight; }
  inline REAL get_weight(){ return weight_; }

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

  virtual int mult_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
  virtual int mult_MH(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
  virtual int mult_MH_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
  
  void* operator new (size_t bytes) { return ::new char[bytes]; }
  void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
  void * operator new(size_t s, void * p) { return p; }

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
  REAL weight_;
  int old_device_;

};

#endif //CUCGMATRIXOPERATOR_H
