#ifndef CUNDARRAY_H
#define CUNDARRAY_H

#include "NDArray.h"
#include "hoNDArray.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

template <class T> class cuNDArray;
template <class T> int cuNDArray_permute(cuNDArray<T>* in, 
					 cuNDArray<T>* out,
					 std::vector<unsigned int> order);

template <class T> class cuNDArray : public NDArray<T>
{
  
  friend int cuNDArray_permute<>(cuNDArray<T>* in, 
				 cuNDArray<T>* out,
				 std::vector<unsigned int> order);
  
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
    
    cudaPointerAttributes attrib;
    if (cudaPointerGetAttributes(&attrib, data) != cudaSuccess) {
      std::cerr << "cuNDArray::create: Unable to determine attributes of pointer" << std::endl;
      return 0;
    }
    
    this->dimensions_ = dimensions;
    this->data_ = data;
    this->device_ = attrib.device;
    this->delete_data_on_destruct_ = delete_data_on_destruct;
    this->elements_ = 1;
    for (unsigned int i = 0; i < this->dimensions_.size(); i++) {
      this->elements_ *= this->dimensions_[i];
    }

    return this->data_;
  }


  virtual int get_device() {return device_;}

  //Copy constructor
  cuNDArray(const cuNDArray<T>& a) {
    this->data_ = 0;
    this->dimensions_ = a.dimensions_;
    if (allocate_memory() == 0) {
      if (a.device_ == this->device_) {
	if (cudaMemcpy(this->data_, a.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice) !=
	    cudaSuccess) {
	  std::cerr << "cuNDArray: Unable to copy data in copy constructor" << std::endl;
	}
      } else {
	//This memory is on a different device, we must move it.
	hoNDArray<T> tmp;
	if (!tmp.create(a.dimensions_)) {
	  std::cerr << "cuNDArray: Unable to allocate temporary memory in copy constructor" << std::endl;
	  deallocate_memory();
	  this->dimensions_.clear();
	}
	
	if (cudaMemcpy(tmp->get_data_ptr(), a.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToHost) !=
	    cudaSuccess) {
	  deallocate_memory();
	  this->data_ = 0;
	  this->dimensions_.clear();
	}

	if (cudaMemcpy(this->data, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) !=
	    cudaSuccess) {
	  deallocate_memory();
	  this->data_ = 0;
	  this->dimensions_.clear();
	}
      }
    }
  }
  
  cuNDArray(hoNDArray<T>& a) {
    this->data_ = 0;
    this->dimensions_ = a.get_dimensions();
    if (allocate_memory() == 0) {
      if (cudaMemcpy(this->data_, a.get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) !=
	  cudaSuccess) {
	deallocate_memory();
	this->data_ = 0;
	this->dimensions_.clear();
      }
    }
  }
  
  
  hoNDArray<T> to_host() {
    hoNDArray<T> ret;
    if (ret.create(this->dimensions_)) {
      if (cudaMemcpy(ret.get_data_ptr(), this->data_, this->elements_*sizeof(T), cudaMemcpyDeviceToHost) !=
	  cudaSuccess) {
	std::cerr << "cuNDArray::to_host(): failed to copy memory from device" << std::endl;
      }
    }
    return ret;
  }


  virtual int permute(std::vector<unsigned int>& dim_order, NDArray<T>* out = 0)
  {

    cuNDArray<T>* out_int = 0;

    //Check ordering array
    if (dim_order.size() > this->dimensions_.size()) {
      std::cerr << "hoNDArray::permute - Invalid length of dimension ordering array" << std::endl;
      return -1;
    }

    std::vector<unsigned int> dim_count(this->dimensions_.size(),0);
    for (unsigned int i = 0; i < dim_order.size(); i++) {
      if (dim_order[i] >= this->dimensions_.size()) {
	std::cerr << "hoNDArray::permute - Invalid dimension order array" << std::endl;
	return -1;
      }
      dim_count[dim_order[i]]++;
    }

    //Create an internal array to store the dimensions
    std::vector<unsigned int> dim_order_int;

    //Check that there are no duplicate dimensions
    for (unsigned int i = 0; i < dim_order.size(); i++) {
      if (dim_count[dim_order[i]] != 1) {
	std::cerr << "hoNDArray::permute - Invalid dimension order array (duplicates)" << std::endl;
	return -1;
      }
      dim_order_int.push_back(dim_order[i]);
    }

    //Pad dimension order array with dimension not mentioned in order array
    if (dim_order_int.size() < this->dimensions_.size()) {
      for (unsigned int i = 0; i < dim_count.size(); i++) {
	if (dim_count[i] == 0) {
	  dim_order_int.push_back(i);
	}
      }
    }

    if (out) {
      out_int = dynamic_cast< cuNDArray<T>* >(out);
      if (!out_int) {
	std::cerr << "cuNDArray::permute: failed to dynamic cast out array pointer" << std::endl;
	return -1;
      }
      for (unsigned int i = 0; i < dim_order_int.size(); i++) {
	if (this->dimensions_[dim_order_int[i]] != out_int->get_size(i)) {
	  std::cerr << "cuNDArray::permute: Dimensions of output array do not match the input array" << std::endl;
	  return -1;
	}
      }
    }

    return cuNDArray_permute(this, out_int, dim_order_int);
  }

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
