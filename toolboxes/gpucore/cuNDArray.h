#ifndef CUNDARRAY_H
#define CUNDARRAY_H
#pragma once

#include "NDArray.h"
#include "hoNDArray.h"

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <boost/shared_ptr.hpp>

template <class T> class cuNDArray;
template <class T> int cuNDArray_permute(cuNDArray<T>* in, cuNDArray<T>* out, std::vector<unsigned int> *order, int shift_mode);

template <class T> class cuNDArray : public NDArray<T>
{
  friend int cuNDArray_permute<>(cuNDArray<T>* in, cuNDArray<T>* out, std::vector<unsigned int> *order, int shift_mode);
  
public:
  
  cuNDArray () : NDArray<T>::NDArray() { cudaGetDevice(&this->device_); }
  
  virtual ~cuNDArray() 
  {
    if (this->delete_data_on_destruct_) {
      deallocate_memory();
    }
  }
  
  // Copy constructor
  cuNDArray(const cuNDArray<T>& a);

  // Constructor from hoNDArray
  cuNDArray(hoNDArray<T> *a);

  // Assignment operator
  cuNDArray& operator=(const cuNDArray<T>& rhs);
  
  virtual T* create(std::vector<unsigned int> *dimensions)
  {
    return NDArray<T>::create(dimensions);
  }

  virtual T* create(std::vector<unsigned int> *dimensions, int device_no)
  {
    if (device_no < 0){
      std::cerr << "cuNDArray::create: illegal device no" << std::endl;
      return 0x0;
    }
    
    int device_no_old;
    if (cudaGetDevice(&device_no_old) != cudaSuccess) {
      std::cerr << "cuNDArray::create: unable to get device no" << std::endl;
      return 0x0;
    }
    
    if (device_no != device_no_old) {
      if (cudaSetDevice(device_no) != cudaSuccess) {
	std::cerr << "cuNDArray::create: unable to set device no" << std::endl;
	return 0x0;
      }
    }

    this->device_ = device_no;
    
    T* ret_value = NDArray<T>::create(dimensions);

    if (device_no != device_no_old) {
      if (cudaSetDevice(device_no_old) != cudaSuccess) {
	std::cerr << "cuNDArray::create unable to restore device no" << std::endl;
      }
    }
    
    return ret_value;
  }
  
  virtual T* create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct = false) 
  {
    if (!data) {
      std::cerr << "cuNDArray::create: 0x0 pointer provided" << std::endl;
      return 0x0;
    }
    
    cudaDeviceProp deviceProp;  
    cudaGetDeviceProperties( &deviceProp, 0);

    if (deviceProp.major >= 2) {
      cudaPointerAttributes attrib;
      if (cudaPointerGetAttributes(&attrib, data) != cudaSuccess) {
	std::cerr << "cuNDArray::create: Unable to determine attributes of pointer" << std::endl;
	return 0x0;
      }
      this->device_ = attrib.device;
    } else {
      cudaGetDevice(&this->device_);
    }
    
    return NDArray<T>::create(dimensions, data, delete_data_on_destruct);
  }

  static boost::shared_ptr< cuNDArray<T> > allocate(std::vector<unsigned int> *dimensions, int device_no) 
  {
    boost::shared_ptr< cuNDArray<T> > ret( new cuNDArray<T> );
    
    if( ret->create(dimensions, device_no) == 0x0 ) {
      std::cerr << "cuNDArray<T>::allocate failed to create array on device " << device_no << std::endl;
      return boost::shared_ptr< cuNDArray<T> >();
    }
    
    return ret;
  }

  static boost::shared_ptr< cuNDArray<T> > allocate(std::vector<unsigned int> *dimensions)
  {
    int tmp_device; cudaGetDevice(&tmp_device);
    return allocate( dimensions, tmp_device );
  }
  
  boost::shared_ptr< hoNDArray<T> > to_host() const {
    boost::shared_ptr< hoNDArray<T> > ret = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>);
    if (ret->create(this->dimensions_.get())) {
      if (cudaMemcpy(ret->get_data_ptr(), this->data_, this->elements_*sizeof(T), cudaMemcpyDeviceToHost) != cudaSuccess) {
	std::cerr << "cuNDArray::to_host(): failed to copy memory from device" << std::endl;
      }
    }
    else{
      std::cerr << "cuNDArray::to_host(): failed to create host array" << std::endl;
    }
    return ret;
  }
  
  virtual int permute(std::vector<unsigned int> *dim_order, NDArray<T> *out = 0, int shift_mode = 0);
  
  virtual int get_device() { return device_; }
  
protected:
  
  int device_;
  
  virtual int allocate_memory()
  {
    deallocate_memory();
    
    this->elements_ = 1;
    for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
      this->elements_ *= (*this->dimensions_)[i];
    } 
    
    size_t size = this->elements_ * sizeof(T);

    if (cudaMalloc((void**) &this->data_,size) != cudaSuccess) {
      std::cerr << "cuNDArray::allocate_memory() : Error allocating CUDA memory" << std::endl;
      this->data_ = 0;
      return -1;
    }
    
    cudaGetDevice(&device_);
    
    return 0;
  }
  
  virtual int deallocate_memory() 
  {
    if (this->data_) {
      if (cudaFree(this->data_) != cudaSuccess) {
	std::cerr << "cuNDArray::deallocate_memory(): failed to delete device memory" << std::endl;
	return -1;
      }
      this->data_ = 0;
    }
    
    return 0;
  }
};

#endif //CUNDARRAY_H
