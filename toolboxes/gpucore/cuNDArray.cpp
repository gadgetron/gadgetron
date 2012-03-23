#include "cuNDArray.h"
#include "vector_td.h"

template <class T> 
cuNDArray<T>::cuNDArray() : NDArray<T>::NDArray() 
{ 
  cudaGetDevice(&this->device_); 
}

template <class T> 
cuNDArray<T>::cuNDArray(hoNDArray<T>  *a) : NDArray<T>::NDArray() 
{
  cudaGetDevice(&this->device_);
  this->dimensions_ = a->get_dimensions();

  if (allocate_memory() == 0) {
    if (cudaMemcpy(this->data_, a->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
      deallocate_memory();
      this->data_ = 0;
      this->dimensions_->clear();
    }
  }
}

template <class T> 
cuNDArray<T>::cuNDArray(const cuNDArray<T>& a) 
{
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

template <class T> 
cuNDArray<T>& cuNDArray<T>::operator=(const cuNDArray<T>& rhs) 
{
  int cur_device; 
  if( cudaGetDevice(&cur_device) != cudaSuccess) {
    std::cerr << "cuNDArray::operator=: unable to get device no" << std::endl;
    return *this;
  }
  
  bool dimensions_match = this->dimensions_equal(&rhs);
  
  if (dimensions_match && (rhs.device_ == cur_device) && (cur_device == this->device_)) {
    
    if (cudaMemcpy(this->data_, rhs.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice) !=
	cudaSuccess) {
      std::cerr << "cuNDArray& operator= failed to copy data" << std::endl;
      return *this;
    }
  } 
  else {
    
    if (cudaSetDevice(this->device_) != cudaSuccess) {
      std::cerr << "cuNDArray::operator=: unable to set device no" << std::endl;
      return *this;
    }
    
    if( !dimensions_match ){
      
      deallocate_memory();
      
      this->elements_ = rhs.elements_;
      this->dimensions_ = rhs.dimensions_;
      
      if (allocate_memory() != 0) {
	std::cerr << "cuNDArray::operator=: failed to allocate memory" << std::endl;
	cudaSetDevice(cur_device);
	return *this;
      }
    }
    
    if (this->device_ == rhs.device_) {
      if (cudaMemcpy(this->data_, rhs.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice) !=
	  cudaSuccess) {
	std::cerr << "cuNDArray::operator=: failed to copy data (2)" << std::endl;
	cudaSetDevice(cur_device);
	return *this;
      }
    } else {

      if( cudaSetDevice(rhs.device_) != cudaSuccess) {
	std::cerr << "cuNDArray::operator=: unable to set device no (2)" << std::endl;
	cudaSetDevice(cur_device);
	return *this;
      }

      boost::shared_ptr< hoNDArray<T> > tmp = rhs.to_host();
      if( cudaSetDevice(this->device_) != cudaSuccess) {
	std::cerr << "cuNDArray::operator=: unable to set device no (3)" << std::endl;
	cudaSetDevice(cur_device);
	return *this;
      }
      
      if (cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
	std::cerr << "cuNDArray::operator=: failed to copy data (3)" << std::endl;
	cudaSetDevice(cur_device);
	return *this;
      }
    }

    if( cudaSetDevice(cur_device) != cudaSuccess) {
      std::cerr << "cuNDArray::operator=: unable to restore to current device" << std::endl;
      return *this;
    }
  }
  
  return *this;
}

template <class T>
cuNDArray<T>::~cuNDArray() 
{ 
  if (this->delete_data_on_destruct_) 
    deallocate_memory();  
}

template <class T>
T* cuNDArray<T>::create(std::vector<unsigned int> *dimensions)
{
  return NDArray<T>::create(dimensions);
}

template <class T>
T* cuNDArray<T>::create(std::vector<unsigned int> *dimensions, int device_no)
{
  if (device_no < 0){
    std::cerr << "cuNDArray::create: illegal device no" << std::endl;
    return 0x0;
  }
  
  this->device_ = device_no; 
  return NDArray<T>::create(dimensions);
}

template <class T>
T* cuNDArray<T>::create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct) 
{
  if (!data) {
    std::cerr << "cuNDArray::create: 0x0 pointer provided" << std::endl;
    return 0x0;
  }
  
  int tmp_device; 
  if( cudaGetDevice(&tmp_device) != cudaSuccess) {
    std::cerr << "cuNDArray::create: Unable to query for device" << std::endl;
    return 0x0;
  }
  
  cudaDeviceProp deviceProp; 
  if( cudaGetDeviceProperties( &deviceProp, tmp_device) != cudaSuccess) {
    std::cerr << "cuNDArray::create: Unable to query device properties" << std::endl;
    return 0x0;
  }
  
  if (deviceProp.unifiedAddressing) {
    cudaPointerAttributes attrib;
    if (cudaPointerGetAttributes(&attrib, data) != cudaSuccess) {
      std::cerr << "cuNDArray::create: Unable to determine attributes of pointer" << std::endl;
      return 0x0;
    }
    this->device_ = attrib.device;
  } else {
    this->device_ = tmp_device;
  }
  
  return NDArray<T>::create(dimensions, data, delete_data_on_destruct);
}

template <class T>
boost::shared_ptr< cuNDArray<T> > cuNDArray<T>::allocate(std::vector<unsigned int> *dimensions)
{
  int tmp_device; cudaGetDevice(&tmp_device);
  return allocate( dimensions, tmp_device );
}

template <class T>
boost::shared_ptr< cuNDArray<T> > cuNDArray<T>::allocate(std::vector<unsigned int> *dimensions, int device_no) 
{
  boost::shared_ptr< cuNDArray<T> > ret( new cuNDArray<T> );
  
  if( ret->create(dimensions, device_no) == 0x0 ) {
    std::cerr << "cuNDArray<T>::allocate failed to create array on device " << device_no << std::endl;
    return boost::shared_ptr< cuNDArray<T> >();
  }
  
  return ret;
}

template <class T>
boost::shared_ptr< hoNDArray<T> > cuNDArray<T>::to_host() const 
{
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

template <class T> 
int cuNDArray<T>::permute(std::vector<unsigned int> *dim_order, NDArray<T> *out, int shift_mode)
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

template <class T> 
int cuNDArray<T>::set_device(int device)
{
  if( device_ == device )
    return 0;

  int cur_device; 
  if( cudaGetDevice(&cur_device) != cudaSuccess) {
    std::cerr << "cuNDArray::set_device: unable to get device no" << std::endl;
    return -1;
  }

  if( cur_device != device_ && cudaSetDevice(device_) != cudaSuccess) {
    std::cerr << "cuNDArray::set_device: unable to set device no" << std::endl;
    return -1;
  }
  
  boost::shared_ptr< hoNDArray<T> > tmp = to_host();

  if( deallocate_memory() != 0 ){
    std::cerr << "cuNDArray::set_device: unable to deallocate memory" << std::endl;
    cudaSetDevice(cur_device);
    return -1;
  }
  
  if( cudaSetDevice(device) != cudaSuccess) {
    std::cerr << "cuNDArray::set_device: unable to set device no (2)" << std::endl;
    cudaSetDevice(cur_device);
    return -1;
  }
 
  device_ = device;
  
  if( allocate_memory() != 0 ) {
    std::cerr << "cuNDArray::set_device: unable to allocate memory" << std::endl;
    cudaSetDevice(cur_device);
    return -1;
  }

  if (cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
    std::cerr << "cuNDArray::set_device: failed to copy data" << std::endl;
    cudaSetDevice(cur_device);
    return -1;
  }

  if( cudaSetDevice(cur_device) != cudaSuccess) {
    std::cerr << "cuNDArray::set_device: unable to restore device to current device" << std::endl;
    return -1;
  }

  return 0;
}

template <class T>
int cuNDArray<T>::allocate_memory()
{
  deallocate_memory();
  
  this->elements_ = 1;
  for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
    this->elements_ *= (*this->dimensions_)[i];
  } 
  
  size_t size = this->elements_ * sizeof(T);
  
  int device_no_old;
  if (cudaGetDevice(&device_no_old) != cudaSuccess) {
    std::cerr << "cuNDArray::allocate_memory: unable to get device no" << std::endl;
    return -1;
  }
  
  if (device_ != device_no_old) {
    if (cudaSetDevice(device_) != cudaSuccess) {
      std::cerr << "cuNDArray::allocate_memory: unable to set device no" << std::endl;
      return -1;
    }
  }
  
  if (cudaMalloc((void**) &this->data_,size) != cudaSuccess) {
    std::cerr << "cuNDArray::allocate_memory() : Error allocating CUDA memory" << std::endl;
    size_t free = 0, total = 0;
    cudaMemGetInfo(&free, &total);
    std::cerr << "CUDA Memory: " << free << " (" << total << ")" << std::endl;
    std::cerr << "   memory requested: " << size << "( " << std::endl;
    for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
      std::cerr << (*this->dimensions_)[i] << " ";
    } 
    std::cerr << " )" << std::endl;
    this->data_ = 0;
    return -1;
  }
  
  if (device_ != device_no_old) {
    if (cudaSetDevice(device_no_old) != cudaSuccess) {
      std::cerr << "cuNDArray::allocate_memory: unable to restore device no" << std::endl;
    }
  }
  
  return 0;
}

template <class T>
int cuNDArray<T>::deallocate_memory()
{
  if (this->data_) {
    
    int device_no_old;
    if (cudaGetDevice(&device_no_old) != cudaSuccess) {
      std::cerr << "cuNDArray::deallocate_memory: unable to get device no" << std::endl;
      return -1;
    }
    
    if (device_ != device_no_old) {
      if (cudaSetDevice(device_) != cudaSuccess) {
	std::cerr << "cuNDArray::deallocate_memory: unable to set device no" << std::endl;
	return -1;
      }
    }
    
    if (cudaFree(this->data_) != cudaSuccess) {
      std::cerr << "cuNDArray::deallocate_memory(): failed to delete device memory" << std::endl;
      return -1;
    }
    
    if (device_ != device_no_old) {
      if (cudaSetDevice(device_no_old) != cudaSuccess) {
	std::cerr << "cuNDArray::allocate_memory: unable to restore device no" << std::endl;
      }
    }
    
    this->data_ = 0;
  }
  
  return 0;
}

//
// Instantiation
//

template EXPORTGPUCORE class cuNDArray< int >;
template EXPORTGPUCORE class cuNDArray< int2 >;
template EXPORTGPUCORE class cuNDArray< int3 >;
template EXPORTGPUCORE class cuNDArray< int4 >;

template EXPORTGPUCORE class cuNDArray< unsigned int >;
template EXPORTGPUCORE class cuNDArray< uint2 >;
template EXPORTGPUCORE class cuNDArray< uint3 >;
template EXPORTGPUCORE class cuNDArray< uint4 >;

template EXPORTGPUCORE class cuNDArray< float >;
template EXPORTGPUCORE class cuNDArray< float2 >;
template EXPORTGPUCORE class cuNDArray< float3 >;
template EXPORTGPUCORE class cuNDArray< float4 >;

template EXPORTGPUCORE class cuNDArray< double >;
template EXPORTGPUCORE class cuNDArray< double2 >;
template EXPORTGPUCORE class cuNDArray< double3 >;
template EXPORTGPUCORE class cuNDArray< double4 >;

template EXPORTGPUCORE class cuNDArray< intd<1>::Type >;
template EXPORTGPUCORE class cuNDArray< intd<2>::Type >;
template EXPORTGPUCORE class cuNDArray< intd<3>::Type >;
template EXPORTGPUCORE class cuNDArray< intd<4>::Type >;

template EXPORTGPUCORE class cuNDArray< uintd<1>::Type >;
template EXPORTGPUCORE class cuNDArray< uintd<2>::Type >;
template EXPORTGPUCORE class cuNDArray< uintd<3>::Type >;
template EXPORTGPUCORE class cuNDArray< uintd<4>::Type >;

template EXPORTGPUCORE class cuNDArray< floatd<1>::Type >;
template EXPORTGPUCORE class cuNDArray< floatd<2>::Type >;
template EXPORTGPUCORE class cuNDArray< floatd<3>::Type >;
template EXPORTGPUCORE class cuNDArray< floatd<4>::Type >;

template EXPORTGPUCORE class cuNDArray< doubled<1>::Type >;
template EXPORTGPUCORE class cuNDArray< doubled<2>::Type >;
template EXPORTGPUCORE class cuNDArray< doubled<3>::Type >;
template EXPORTGPUCORE class cuNDArray< doubled<4>::Type >;

template EXPORTGPUCORE class cuNDArray< intd1 >;
template EXPORTGPUCORE class cuNDArray< intd2 >;
template EXPORTGPUCORE class cuNDArray< intd3 >;
template EXPORTGPUCORE class cuNDArray< intd4 >;

template EXPORTGPUCORE class cuNDArray< uintd1 >;
template EXPORTGPUCORE class cuNDArray< uintd2 >;
template EXPORTGPUCORE class cuNDArray< uintd3 >;
template EXPORTGPUCORE class cuNDArray< uintd4 >;

template EXPORTGPUCORE class cuNDArray< floatd1 >;
template EXPORTGPUCORE class cuNDArray< floatd2 >;
template EXPORTGPUCORE class cuNDArray< floatd3 >;
template EXPORTGPUCORE class cuNDArray< floatd4 >;

template EXPORTGPUCORE class cuNDArray< doubled1 >;
template EXPORTGPUCORE class cuNDArray< doubled2 >;
template EXPORTGPUCORE class cuNDArray< doubled3 >;
template EXPORTGPUCORE class cuNDArray< doubled4 >;

template EXPORTGPUCORE class cuNDArray<float_complext>;
template EXPORTGPUCORE class cuNDArray<double_complext>;
