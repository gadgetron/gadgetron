#pragma once
#include "cuCGMatrixOperator.h"

template <class REAL, class T> class cuCGMatrixOperatorDevice : public cuCGMatrixOperator<REAL,T>
{
public:
  cuCGMatrixOperatorDevice( int device = -1 ) : cuCGMatrixOperator<REAL,T>() 
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
  
  virtual ~cuCGMatrixOperatorDevice() {}
  
  virtual int mult_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false){
    int ret = set_device();
    if( ret == 0 ) ret = mult_M_device( in, out, accumulate );
    if( ret == 0 ) ret = restore_device();
    return ret;
  }

  virtual int mult_MH(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false){
    int ret = set_device();
    if( ret == 0 ) ret = mult_MH_device( in, out, accumulate );
    if( ret == 0 ) ret = restore_device();
    return ret;
  }

  virtual int mult_MH_M(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false){
    int ret = set_device();
    if( ret == 0 ) ret = mult_MH_M_device( in, out, accumulate );
    if( ret == 0 ) ret = restore_device();
    return ret;
  }
  
  virtual int mult_M_device(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
  virtual int mult_MH_device(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;
  virtual int mult_MH_M_device(cuNDArray<T>* in, cuNDArray<T>* out, bool accumulate = false) = 0;

  void* operator new (size_t bytes) { return ::new char[bytes]; }
  void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
  void * operator new(size_t s, void * p) { return p; }

protected:
  virtual int set_device()
  {
    if( cudaGetDevice( &old_device_ ) != cudaSuccess ){
      std::cerr << "cuCGMatrixOperatorDevice::set_device: unable to get current device." << std::endl ;
      return -1;
    }
    if( device_ != old_device_ && cudaSetDevice(device_) != cudaSuccess) {
      std::cerr << "cuCGMatrixOperatorDevice::set_device: unable to set device " << device_ << std::endl;
      return -1;
    }
    return 0;
  }

  virtual int restore_device()
  {
    if( device_ != old_device_ && cudaSetDevice(old_device_) != cudaSuccess) {
      std::cerr << "cuCGMatrixOperatorDevice::restore_device: unable to set device " << old_device_ << std::endl;
      return -1;
    }
    return 0;
  }

protected:
  int device_;

private:
  int old_device_;
};
