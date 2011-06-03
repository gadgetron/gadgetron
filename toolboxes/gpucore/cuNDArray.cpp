#include "cuNDArray.h"
#include "vector_td.h"

template <class T> cuNDArray<T>::cuNDArray(hoNDArray<T>  *a) {
  this->data_ = 0;
  this->dimensions_ = a->get_dimensions();
  if (allocate_memory() == 0) {
    if (cudaMemcpy(this->data_, a->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) !=
	cudaSuccess) {
      deallocate_memory();
      this->data_ = 0;
      this->dimensions_->clear();
    }
  }
}

template <class T> cuNDArray<T>::cuNDArray(const cuNDArray<T>& a) {
  cudaGetDevice(&this->device_);
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
      cudaSetDevice(a.device_);
      boost::shared_ptr< hoNDArray<T> > tmp = a.to_host();
      cudaSetDevice(this->device_);
      
      if (cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) !=cudaSuccess) {
	deallocate_memory();
	this->data_ = 0;
	this->dimensions_->clear();
      }
    }
  }
}

template <class T> cuNDArray<T>&  cuNDArray<T>::operator=(const cuNDArray<T>& rhs) {
  int old_device = this->device_;
  int cur_device;
  cudaGetDevice(&cur_device);
  
  bool dimensions_match = this->dimensions_equal(&rhs);
  
  if (dimensions_match &&
      (rhs.device_ == cur_device) &&
      (cur_device == old_device)) {

    if (cudaMemcpy(this->data_, rhs.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice) !=
	cudaSuccess) {
      std::cerr << "cuNDArray& operator= failed to copy data" << std::endl;
      return *this;
	}
  } else {

    if (old_device != cur_device) {
      cudaSetDevice(old_device);
    }
    deallocate_memory();
    if (old_device != cur_device) {
      cudaSetDevice(cur_device);
    }
    
    this->elements_ = rhs.elements_;
    this->dimensions_ = rhs.dimensions_;
    if (allocate_memory() != 0) {
      std::cerr << "cuNDArray& operator= failed to allocate memory" << std::endl;
      return *this;
    }
    
    if (this->device_ == rhs.device_) {
      if (cudaMemcpy(this->data_, rhs.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice) !=
	  cudaSuccess) {
	std::cerr << "cuNDArray& operator= failed to copy data (2)" << std::endl;
	return *this;
      }
    } else {
      cudaSetDevice(rhs.device_);
      boost::shared_ptr< hoNDArray<T> > tmp = rhs.to_host();
      cudaSetDevice(this->device_);
      
      if (cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
	std::cerr << "cuNDArray& operator= failed to copy data (3)" << std::endl;
	return *this;
      }
    }
  }
  
  return *this;
}

template <class T> int cuNDArray<T>::permute(std::vector<unsigned int> *dim_order, NDArray<T> *out, int shift_mode)
{

  cuNDArray<T>* out_int = 0x0;
  
  //Check ordering array
  if (dim_order->size() > this->dimensions_->size()) {
    std::cerr << "hoNDArray::permute - Invalid length of dimension ordering array" << std::endl;
    return -1;
  }
  
  std::vector<unsigned int> dim_count(this->dimensions_->size(),0);
  for (unsigned int i = 0; i < dim_order->size(); i++) {
    if ((*dim_order)[i] >= this->dimensions_->size()) {
      std::cerr << "hoNDArray::permute - Invalid dimension order array" << std::endl;
      return -1;
    }
    dim_count[(*dim_order)[i]]++;
  }
  
  //Create an internal array to store the dimensions
  std::vector<unsigned int> dim_order_int;
  
  //Check that there are no duplicate dimensions
  for (unsigned int i = 0; i < dim_order->size(); i++) {
    if (dim_count[(*dim_order)[i]] != 1) {
      std::cerr << "hoNDArray::permute - Invalid dimension order array (duplicates)" << std::endl;
      return -1;
    }
    dim_order_int.push_back((*dim_order)[i]);
  }
  
  //Pad dimension order array with dimension not mentioned in order array
  if (dim_order_int.size() < this->dimensions_->size()) {
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
      if ((*this->dimensions_)[dim_order_int[i]] != out_int->get_size(i)) {
	std::cerr << "cuNDArray::permute: Dimensions of output array do not match the input array" << std::endl;
	return -1;
      }
    }
  }
  
  return cuNDArray_permute(this, out_int, &dim_order_int, shift_mode);
}

//
// Instantiation
//

template class cuNDArray< int >;
template class cuNDArray< int2 >;
template class cuNDArray< int3 >;
template class cuNDArray< int4 >;

template class cuNDArray< unsigned int >;
template class cuNDArray< uint2 >;
template class cuNDArray< uint3 >;
template class cuNDArray< uint4 >;

template class cuNDArray< float >;
template class cuNDArray< float2 >;
template class cuNDArray< float3 >;
template class cuNDArray< float4 >;

template class cuNDArray< double >;
template class cuNDArray< double2 >;
template class cuNDArray< double3 >;
template class cuNDArray< double4 >;

template class cuNDArray< intd<1>::Type >;
template class cuNDArray< intd<2>::Type >;
template class cuNDArray< intd<3>::Type >;
template class cuNDArray< intd<4>::Type >;

template class cuNDArray< uintd<1>::Type >;
template class cuNDArray< uintd<2>::Type >;
template class cuNDArray< uintd<3>::Type >;
template class cuNDArray< uintd<4>::Type >;

template class cuNDArray< floatd<1>::Type >;
template class cuNDArray< floatd<2>::Type >;
template class cuNDArray< floatd<3>::Type >;
template class cuNDArray< floatd<4>::Type >;

template class cuNDArray< doubled<1>::Type >;
template class cuNDArray< doubled<2>::Type >;
template class cuNDArray< doubled<3>::Type >;
template class cuNDArray< doubled<4>::Type >;

template class cuNDArray< intd1 >;
template class cuNDArray< intd2 >;
template class cuNDArray< intd3 >;
template class cuNDArray< intd4 >;

template class cuNDArray< uintd1 >;
template class cuNDArray< uintd2 >;
template class cuNDArray< uintd3 >;
template class cuNDArray< uintd4 >;

template class cuNDArray< floatd1 >;
template class cuNDArray< floatd2 >;
template class cuNDArray< floatd3 >;
template class cuNDArray< floatd4 >;

template class cuNDArray< doubled1 >;
template class cuNDArray< doubled2 >;
template class cuNDArray< doubled3 >;
template class cuNDArray< doubled4 >;

template class cuNDArray<float_complext>;
template class cuNDArray<double_complext>;
