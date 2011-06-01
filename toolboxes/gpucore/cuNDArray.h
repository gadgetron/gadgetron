#ifndef CUNDARRAY_H
#define CUNDARRAY_H

#pragma once
#include "gadgetron_export.h"

#include "NDArray.h"
#include "hoNDArray.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

template <class T> class EXPORTGPUCORE cuNDArray;
template <class T> EXPORTGPUCORE int cuNDArray_permute(cuNDArray<T>* in, cuNDArray<T>* out, std::vector<unsigned int> order, int shift_mode);

template <class T> class EXPORTGPUCORE cuNDArray : public NDArray<T>
{
  
  friend EXPORTGPUCORE int cuNDArray_permute<>(cuNDArray<T>* in, 
				 cuNDArray<T>* out,
				 std::vector<unsigned int> order,
				 int shift_mode);
  
 public:
  cuNDArray () 
    : NDArray<T>::NDArray()
    , device_(0)
  {
    
  }

  ~cuNDArray() {
    if (this->delete_data_on_destruct_) {
      deallocate_memory();
    }
  }

  static cuNDArray<T>* allocate(std::vector<unsigned int> dimensions) {
    cuNDArray<T>* ret = 0;
    
    ret = new cuNDArray<T>;

    if (ret) {
      if (!ret->create(dimensions)) {
	std::cerr << "cuNDArray<T>* allocate failed to allocate memory in array" << std::endl;
	delete ret;
	ret = 0;
      }
    }
    
    return ret;
  }

  static cuNDArray<T>* allocate(std::vector<unsigned int> dimensions, int device_no) {
    cuNDArray<T>* ret = 0;
    
    ret = new cuNDArray<T>;

    if (ret) {
      if (!ret->create(dimensions, device_no)) {
	std::cerr << "cuNDArray<T>* allocate failed to allocate memory in array" << std::endl;
	delete ret;
	ret = 0;
      }
    }
    
    return ret;
  }


  virtual T* create(std::vector<unsigned int> dimensions) {
    this->dimensions_ = dimensions; 
    allocate_memory();
    return this->get_data_ptr();
  }


  virtual T* create(std::vector<unsigned int>& dimensions, int device_no) {
    int device_no_old;
    if (cudaGetDevice(&device_no_old) != cudaSuccess) {
      std::cerr << "cuNDArray::create unable to get device no" << std::endl;
    }

    if ((device_no >= 0) && (device_no != device_no_old) ) {
      if (cudaSetDevice(device_no) != cudaSuccess) {
	std::cerr << "cuNDArray::create unable to set device no" << std::endl;
      }
    }

    T* ret_value = create(dimensions);

    if ((device_no >= 0) && (device_no != device_no_old) ) {
      if (cudaSetDevice(device_no_old) != cudaSuccess) {
	std::cerr << "cuNDArray::create unable to set device no" << std::endl;
      }
    }

    return ret_value;
  }


  virtual T* create(std::vector<unsigned int> dimensions, T* data, 
		    bool delete_data_on_destruct = false) 
  {
    if (!data) {
      std::cerr << "cuNDArray::create : null pointer provided for data" << std::endl;
      return 0;
    }
    

    cudaDeviceProp deviceProp;  
    cudaGetDeviceProperties( &deviceProp, 0);

    if (deviceProp.major >= 2) {
      cudaPointerAttributes attrib;
      if (cudaPointerGetAttributes(&attrib, data) != cudaSuccess) {
	std::cerr << "cuNDArray::create: Unable to determine attributes of pointer" << std::endl;
	return 0;
      }
      this->device_ = attrib.device;
    } else {
      cudaGetDevice(&this->device_);
    }
    
    this->dimensions_ = dimensions;
    this->data_ = data;
    this->delete_data_on_destruct_ = delete_data_on_destruct;
    this->elements_ = 1;
    for (unsigned int i = 0; i < this->dimensions_.size(); i++) {
      this->elements_ *= this->dimensions_[i];
    }

    return this->data_;
  }


  virtual int get_device() {return device_;}

  //Copy constructor
  cuNDArray(const cuNDArray<T>& a);

  //Assignment operator
  cuNDArray& operator=(const cuNDArray<T>& rhs);

  cuNDArray(hoNDArray<T>& a);

  
  hoNDArray<T> to_host() const {
    hoNDArray<T> ret;
    if (ret.create(this->dimensions_)) {
      if (cudaMemcpy(ret.get_data_ptr(), this->data_, this->elements_*sizeof(T), cudaMemcpyDeviceToHost) !=
	  cudaSuccess) {
	std::cerr << "cuNDArray::to_host(): failed to copy memory from device" << std::endl;
      }
    }
    return ret;
  }


  virtual int permute(std::vector<unsigned int>& dim_order, 
		      NDArray<T>* out = 0,
		      int shift_mode = 0);

 protected:
  int device_;

  virtual int allocate_memory()
  {
    deallocate_memory();
    
    this->elements_ = 1;
    for (unsigned int i = 0; i < this->dimensions_.size(); i++) {
      this->elements_ *= this->dimensions_[i];
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
  
  virtual int deallocate_memory() {
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
