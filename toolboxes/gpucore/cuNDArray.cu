#include "cuNDArray.h"
#include "vector_td.h"
#include <thrust/iterator/constant_iterator.h>
#include <stdexcept>
#include <sstream>
#include "vector_td_utilities.h"
#include "check_CUDA.h"

template <class T> 
cuNDArray<T>::cuNDArray() : NDArray<T>::NDArray() 
{ 
  cudaGetDevice(&this->device_); 
}

template <class T> 
cuNDArray<T>::cuNDArray(std::vector<unsigned int> *dimensions) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions);
}

template <class T>
cuNDArray<T>::cuNDArray(std::vector<unsigned int> *dimensions, int device_no) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions,device_no);
}

template <class T>
cuNDArray<T>::cuNDArray(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions,data,delete_data_on_destruct);
}


template <class T>
cuNDArray<T>::cuNDArray(boost::shared_ptr<std::vector<unsigned int> > dimensions) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions.get());
}

template <class T>
cuNDArray<T>::cuNDArray(boost::shared_ptr<std::vector<unsigned int> > dimensions, int device_no) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions.get(),device_no);
}

template <class T>
cuNDArray<T>::cuNDArray(boost::shared_ptr<std::vector<unsigned int> > dimensions, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions.get(),data,delete_data_on_destruct);
}

template <class T>
cuNDArray<T>::cuNDArray(hoNDArray<T> *a) : NDArray<T>::NDArray() 
{
  cudaGetDevice(&this->device_);
  this->dimensions_ = a->get_dimensions();

  allocate_memory();
	if (cudaMemcpy(this->data_, a->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
	  deallocate_memory();
	  this->data_ = 0;
	  this->dimensions_->clear();

  }
}

template <class T> 
cuNDArray<T>::cuNDArray(const cuNDArray<T>& a) 
{
  cudaGetDevice(&this->device_);
  this->data_ = 0;
  this->dimensions_ = a.dimensions_;

  allocate_memory();
	if (a.device_ == this->device_) {
	  if (cudaMemcpy(this->data_, a.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice) !=
	  cudaSuccess) {
		  throw cuda_error("cuNDArray: Unable to copy data in copy constructor");
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

template <class T> 
cuNDArray<T>& cuNDArray<T>::operator=(const cuNDArray<T>& rhs)
{
  int cur_device; 
  if( cudaGetDevice(&cur_device) != cudaSuccess) {
    throw cuda_error("cuNDArray::operator=: unable to get device no");

  }
  
  bool dimensions_match = this->dimensions_equal(&rhs);
  
  if (dimensions_match && (rhs.device_ == cur_device) && (cur_device == this->device_)) {
    
    if (cudaMemcpy(this->data_, rhs.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice) !=
	cudaSuccess) {
      throw cuda_error("cuNDArray& operator= failed to copy data");

    }
  } 
  else {
    
    if (cudaSetDevice(this->device_) != cudaSuccess) {
      throw cuda_error("cuNDArray::operator=: unable to set device no");

    }
    
    if( !dimensions_match ){
      
      deallocate_memory();
      
      this->elements_ = rhs.elements_;
      this->dimensions_ = rhs.dimensions_;
      
      allocate_memory();

    }
    
    if (this->device_ == rhs.device_) {
      if (cudaMemcpy(this->data_, rhs.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice) !=
	  cudaSuccess) {

	cudaSetDevice(cur_device);
	throw cuda_error("cuNDArray::operator=: failed to copy data (2)");

      }
    } else {

      if( cudaSetDevice(rhs.device_) != cudaSuccess) {

	cudaSetDevice(cur_device);
	throw cuda_error("cuNDArray::operator=: unable to set device no (2)");

      }

      boost::shared_ptr< hoNDArray<T> > tmp = rhs.to_host();
      if( cudaSetDevice(this->device_) != cudaSuccess) {

	cudaSetDevice(cur_device);
	throw cuda_error("cuNDArray::operator=: unable to set device no (3)");

      }
      
      if (cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {

	cudaSetDevice(cur_device);
	throw cuda_error("cuNDArray::operator=: failed to copy data (3)");

      }
    }

    if( cudaSetDevice(cur_device) != cudaSuccess) {
      throw cuda_error("cuNDArray::operator=: unable to restore to current device");

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
    throw cuda_error("cuNDArray::create: illegal device no");

  }
  
  this->device_ = device_no; 
  return NDArray<T>::create(dimensions);
}

template <class T>
T* cuNDArray<T>::create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct) 
{
  if (!data) {
    throw std::logic_error("cuNDArray::create: 0x0 pointer provided");
  }
  
  int tmp_device; 
  if( cudaGetDevice(&tmp_device) != cudaSuccess) {
    throw cuda_error("cuNDArray::create: Unable to query for device");

  }
  
  cudaDeviceProp deviceProp; 
  if( cudaGetDeviceProperties( &deviceProp, tmp_device) != cudaSuccess) {
    throw cuda_error("cuNDArray::create: Unable to query device properties");

  }
  
  if (deviceProp.unifiedAddressing) {
    cudaPointerAttributes attrib;
    if (cudaPointerGetAttributes(&attrib, data) != cudaSuccess) {
      throw cuda_error("cuNDArray::create: Unable to determine attributes of pointer");

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
    throw cuda_error("cuNDArray<T>::allocate failed to create array on device ");

  }
  
  return ret;
}

template <class T>
boost::shared_ptr< hoNDArray<T> > cuNDArray<T>::to_host() const 
{
  boost::shared_ptr< hoNDArray<T> > ret = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>);
  ret->create(this->dimensions_.get());
    if (cudaMemcpy(ret->get_data_ptr(), this->data_, this->elements_*sizeof(T), cudaMemcpyDeviceToHost) != cudaSuccess) {
      throw cuda_error("cuNDArray::to_host(): failed to copy memory from device");
    }

  return ret;
}

template <class T> 
void cuNDArray<T>::permute(std::vector<unsigned int> *dim_order, NDArray<T> *out, int shift_mode)
{

  cuNDArray<T>* out_int = 0x0;
  
  //Check ordering array
  if (dim_order->size() > this->dimensions_->size()) {
    throw std::logic_error("hoNDArray::permute - Invalid length of dimension ordering array");
  }
  
  std::vector<unsigned int> dim_count(this->dimensions_->size(),0);
  for (unsigned int i = 0; i < dim_order->size(); i++) {
    if ((*dim_order)[i] >= this->dimensions_->size()) {
      throw std::logic_error("hoNDArray::permute - Invalid dimension order array");
    }
    dim_count[(*dim_order)[i]]++;
  }
  
  //Create an internal array to store the dimensions
  std::vector<unsigned int> dim_order_int;
  
  //Check that there are no duplicate dimensions
  for (unsigned int i = 0; i < dim_order->size(); i++) {
    if (dim_count[(*dim_order)[i]] != 1) {
      std::logic_error("hoNDArray::permute - Invalid dimension order array (duplicates)");
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
      std::runtime_error("cuNDArray::permute: failed to dynamic cast out array pointer");
    }
    for (unsigned int i = 0; i < dim_order_int.size(); i++) {
      if ((*this->dimensions_)[dim_order_int[i]] != out_int->get_size(i)) {
	std::logic_error("cuNDArray::permute: Dimensions of output array do not match the input array");
      }
    }
  }
  
  cuNDArray_permute(this, out_int, &dim_order_int, shift_mode);
}

template <class T> 
void cuNDArray<T>::set_device(int device)
{
  if( device_ == device )
    return;

  int cur_device; 
  if( cudaGetDevice(&cur_device) != cudaSuccess) {
    throw cuda_error("cuNDArray::set_device: unable to get device no");

  }

  if( cur_device != device_ && cudaSetDevice(device_) != cudaSuccess) {
    throw cuda_error("cuNDArray::set_device: unable to set device no");

  }
  
  boost::shared_ptr< hoNDArray<T> > tmp = to_host();

  deallocate_memory();

  
  if( cudaSetDevice(device) != cudaSuccess) {
	  cudaSetDevice(cur_device);
	  throw cuda_error("cuNDArray::set_device: unable to set device no (2)");
  }
 
  device_ = device;
  
  allocate_memory();
  if (cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
	  cudaSetDevice(cur_device);
    throw cuda_error("cuNDArray::set_device: failed to copy data");
  }

  if( cudaSetDevice(cur_device) != cudaSuccess) {
    throw cuda_error("cuNDArray::set_device: unable to restore device to current device");
  }

}

template <class T>
void cuNDArray<T>::allocate_memory()
{
  deallocate_memory();
  
  this->elements_ = 1;
  if (this->dimensions_->empty())
  	throw std::runtime_error("cuNDArray::allocate_memory() : dimensions is empty.");
  for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
    this->elements_ *= (*this->dimensions_)[i];
  } 
  
  size_t size = this->elements_ * sizeof(T);
  
  int device_no_old;
  if (cudaGetDevice(&device_no_old) != cudaSuccess) {
    throw cuda_error("cuNDArray::allocate_memory: unable to get device no");
  }
  
  if (device_ != device_no_old) {
    if (cudaSetDevice(device_) != cudaSuccess) {
      throw cuda_error("cuNDArray::allocate_memory: unable to set device no");
    }
  }
  
  if (cudaMalloc((void**) &this->data_,size) != cudaSuccess) {
	  size_t free = 0, total = 0;
	  cudaMemGetInfo(&free, &total);
    std::stringstream err("cuNDArray::allocate_memory() : Error allocating CUDA memory");
    err << "CUDA Memory: " << free << " (" << total << ")";

    err << "   memory requested: " << size << "( ";
    for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
      std::cerr << (*this->dimensions_)[i] << " ";
    } 
    err << ")";
    this->data_ = 0;
    throw cuda_error(err.str());

  }
  
  if (device_ != device_no_old) {
    if (cudaSetDevice(device_no_old) != cudaSuccess) {
      throw cuda_error("cuNDArray::allocate_memory: unable to restore device no");
    }
  }

}

template <class T>
void cuNDArray<T>::deallocate_memory()
{
  if (this->data_) {
    
    int device_no_old;
    if (cudaGetDevice(&device_no_old) != cudaSuccess) {
      throw cuda_error("cuNDArray::deallocate_memory: unable to get device no");

    }
    
    if (device_ != device_no_old) {
      if (cudaSetDevice(device_) != cudaSuccess) {
	throw cuda_error("cuNDArray::deallocate_memory: unable to set device no");

      }
    }
    
    if (cudaFree(this->data_) != cudaSuccess) {
      throw cuda_error("cuNDArray::deallocate_memory(): failed to delete device memory");

    }
    
    if (device_ != device_no_old) {
      if (cudaSetDevice(device_no_old) != cudaSuccess) {
    	  throw cuda_error("cuNDArray::allocate_memory: unable to restore device no");
      }
    }
    
    this->data_ = 0;
  }
  

}
template<typename T>
struct cuNDA_abs : public thrust::unary_function<T,T>
{
 __host__ __device__ T operator()(const T &x) const {return abs(x);}
};

template<typename T, unsigned int D>
struct cuNDA_abs<vector_td<T,D> > : public thrust::unary_function<vector_td<T,D>, vector_td<T,D> >
{
 __host__ __device__ vector_td<T,D> operator()(const vector_td<T,D> &x) const {
	 vector_td<T,D> res;
	 for (int i = 0; i < D; i++) res[i] = abs(x[i]);
	 return res;
 	 }
};

template <class T>
void cuNDArray<T>::abs()
{
	thrust::device_ptr<T> devPtr = this->get_device_ptr();
	thrust::transform(devPtr,devPtr+this->get_number_of_elements(),devPtr,cuNDA_abs<T>());
}



template<typename T>
struct cuNDA_sqrt : public thrust::unary_function<T,T>
{
 __host__ __device__ T operator()(const T &x) const {return sqrt(x);}
};

template<typename T, unsigned int D>
struct cuNDA_sqrt<vector_td<T,D> > : public thrust::unary_function<vector_td<T,D>, vector_td<T,D> >
{
 __host__ __device__ vector_td<T,D> operator()(const vector_td<T,D> &x) const {
	 vector_td<T,D> res;
	 for (int i = 0; i < D; i++) res[i] = sqrt(x[i]);
	 return res;
 	 }
};

template <class T>
void cuNDArray<T>::sqrt()
{
	thrust::device_ptr<T> devPtr = this->get_device_ptr();
	thrust::transform(devPtr,devPtr+this->get_number_of_elements(),devPtr,cuNDA_sqrt<T>());
}



template <class T>
void cuNDArray<T>::clear()
{
	cudaMemset(this->get_data_ptr(),0,sizeof(T)*this->get_number_of_elements());
}

template <class T>
void cuNDArray<T>::fill(T val)
{
	thrust::device_ptr<T> devPtr = this->get_device_ptr();
	thrust::fill(devPtr,devPtr+this->get_number_of_elements(),val);
}

template<typename T>
struct cuNDA_reciprocal : public thrust::unary_function<T,T>
{
 __host__ __device__ T operator()(const T &x) const {return 1/x;}
};

template<typename T, unsigned int D>
struct cuNDA_reciprocal<vector_td<T,D> > : public thrust::unary_function<vector_td<T,D>, vector_td<T,D> >
{
 __host__ __device__ vector_td<T,D> operator()(const vector_td<T,D> &x) const {
	 vector_td<T,D> res;
	 for (int i = 0; i < D; i++) res[i] = 1/(x[i]);
	 return res;
 	 }
};
template <class T>
void cuNDArray<T>::reciprocal()
{
	thrust::device_ptr<T> devPtr = this->get_device_ptr();
	thrust::transform(devPtr,devPtr+this->get_number_of_elements(),devPtr,cuNDA_reciprocal<T>());
}



template<class T> void operator+= (cuNDArray<T> & x , cuNDArray<T> & y){
	thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), thrust::plus<T>());
}

template<class T> void operator+= (cuNDArray<T> & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::plus<T>());
}

template<class T> void operator*= (cuNDArray<T> & x , cuNDArray<T> & y){
	thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), thrust::multiplies<T>());
}

template<class T> void operator*= (cuNDArray<T> & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::multiplies<T>());
}

template<class T> void operator-= (cuNDArray<T> & x , cuNDArray<T> & y){
	thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), thrust::minus<T>());
}

template<class T> void operator-= (cuNDArray<T> & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::minus<T>());
}
template<class T> void operator/= (cuNDArray<T> & x , cuNDArray<T> & y){
	thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), thrust::divides<T>());
}

template<class T> void operator/= (cuNDArray<T> & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::divides<T>());
}

template<typename T>
struct cuNDA_plus : public thrust::binary_function<complext<T> ,T,complext<T> >
{
 __host__ __device__ complext<T> operator()(const complext<T> &x, const T &y) const {return x+y;}

};

template<typename T>
struct cuNDA_minus : public thrust::binary_function<complext<T> ,T,complext<T> >
{
 __host__ __device__ complext<T> operator()(const complext<T> &x, const T &y) const {return x-y;}
};

template<typename T>
struct cuNDA_multiplies : public thrust::binary_function<complext<T> ,T,complext<T> >
{
 __host__ __device__ complext<T> operator()(const complext<T> &x, const T &y) const {return x*y;}
};

template<typename T>
struct cuNDA_divides : public thrust::binary_function<complext<T> ,T,complext<T> >
{
 __host__ __device__ complext<T> operator()(const complext<T> &x, const T &y) const {return x/y;}
};

template<class T> void operator+= (cuNDArray< complext<T> > & x , cuNDArray<T> & y){
	thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), cuNDA_plus<T>());
}

template<class T> void operator+= (cuNDArray<complext<T> > & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_plus<T>());
}

template<class T> void operator*= (cuNDArray< complext<T> > & x , cuNDArray<T> & y){
	thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), cuNDA_multiplies<T>());
}

template<class T> void operator*= (cuNDArray<complext<T> > & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_multiplies<T>());
}

template<class T> void operator-= (cuNDArray< complext<T> > & x , cuNDArray<T> & y){
	thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), cuNDA_minus<T>());
}

template<class T> void operator-= (cuNDArray<complext<T> > & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_minus<T>());
}
template<class T> void operator/= (cuNDArray< complext<T> > & x , cuNDArray<T> & y){
	thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), cuNDA_divides<T>());
}

template<class T> void operator/= (cuNDArray<complext<T> > & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_divides<T>());
}


//
// Instantiation
//

template EXPORTGPUCORE class cuNDArray< int >;
/*
template EXPORTGPUCORE class cuNDArray< int2 >;
template EXPORTGPUCORE class cuNDArray< int3 >;
template EXPORTGPUCORE class cuNDArray< int4 >;
*/
template void operator+=<int> (cuNDArray<int> & x , cuNDArray<int> & y);
template void operator-=<int> (cuNDArray<int> & x , cuNDArray<int> & y);
template void operator*=<int> (cuNDArray<int> & x , cuNDArray<int> & y);
template void operator/=<int> (cuNDArray<int> & x , cuNDArray<int> & y);

template void operator+=<int> (cuNDArray<int> & x , int y);
template void operator-=<int> (cuNDArray<int> & x , int y);
template void operator*=<int> (cuNDArray<int> & x , int y);
template void operator/=<int> (cuNDArray<int> & x , int y);

template EXPORTGPUCORE class cuNDArray< unsigned int >;

/*
template EXPORTGPUCORE class cuNDArray< uint2 >;
template EXPORTGPUCORE class cuNDArray< uint3 >;
template EXPORTGPUCORE class cuNDArray< uint4 >;
*/

template void operator+=<unsigned int> (cuNDArray<unsigned int> & x , cuNDArray<unsigned int> & y);
template void operator-=<unsigned int> (cuNDArray<unsigned int> & x , cuNDArray<unsigned int> & y);
template void operator*=<unsigned int> (cuNDArray<unsigned int> & x , cuNDArray<unsigned int> & y);
template void operator/=<unsigned int> (cuNDArray<unsigned int> & x , cuNDArray<unsigned int> & y);

template void operator+=<unsigned int> (cuNDArray<unsigned int> & x , unsigned int y);
template void operator-=<unsigned int> (cuNDArray<unsigned int> & x , unsigned int y);
template void operator*=<unsigned int> (cuNDArray<unsigned int> & x , unsigned int y);
template void operator/=<unsigned int> (cuNDArray<unsigned int> & x , unsigned int y);

template EXPORTGPUCORE class cuNDArray< float >;



/*template EXPORTGPUCORE class cuNDArray< float2 >;
template EXPORTGPUCORE class cuNDArray< float3 >;
template EXPORTGPUCORE class cuNDArray< float4 >;
*/

template void operator+=<float> (cuNDArray<float> & x , cuNDArray<float> & y);
template void operator-=<float> (cuNDArray<float> & x , cuNDArray<float> & y);
template void operator*=<float> (cuNDArray<float> & x , cuNDArray<float> & y);
template void operator/=<float> (cuNDArray<float> & x , cuNDArray<float> & y);

template void operator+=<float> (cuNDArray<float> & x , float y);
template void operator-=<float> (cuNDArray<float> & x , float y);
template void operator*=<float> (cuNDArray<float> & x , float y);
template void operator/=<float> (cuNDArray<float> & x , float y);


template void operator+=<float> (cuNDArray<float_complext> & x , cuNDArray<float> & y);
template void operator-=<float> (cuNDArray<float_complext> & x , cuNDArray<float> & y);
template void operator*=<float> (cuNDArray<float_complext> & x , cuNDArray<float> & y);
template void operator/=<float> (cuNDArray<float_complext> & x , cuNDArray<float> & y);


template void operator+=<float> (cuNDArray<float_complext> & x , float y);
template void operator-=<float> (cuNDArray<float_complext> & x , float y);
template void operator*=<float> (cuNDArray<float_complext> & x , float y);
template void operator/=<float> (cuNDArray<float_complext> & x , float y);


template EXPORTGPUCORE class cuNDArray< double >;

/*
template EXPORTGPUCORE class cuNDArray< double2 >;
template EXPORTGPUCORE class cuNDArray< double3 >;
template EXPORTGPUCORE class cuNDArray< double4 >;
*/

template void operator+=<double> (cuNDArray<double> & x , cuNDArray<double> & y);
template void operator-=<double> (cuNDArray<double> & x , cuNDArray<double> & y);
template void operator*=<double> (cuNDArray<double> & x , cuNDArray<double> & y);
template void operator/=<double> (cuNDArray<double> & x , cuNDArray<double> & y);

template void operator+=<double> (cuNDArray<double> & x , double y);
template void operator-=<double> (cuNDArray<double> & x , double y);
template void operator*=<double> (cuNDArray<double> & x , double y);
template void operator/=<double> (cuNDArray<double> & x , double y);


template EXPORTGPUCORE void operator+=<double> (cuNDArray<double_complext> & x , cuNDArray<double> & y);
template EXPORTGPUCORE void operator-=<double> (cuNDArray<double_complext> & x , cuNDArray<double> & y);
template EXPORTGPUCORE void operator*=<double> (cuNDArray<double_complext> & x , cuNDArray<double> & y);
template EXPORTGPUCORE void operator /=<double> (cuNDArray<double_complext> & x , cuNDArray<double> & y);

template void operator+=<double> (cuNDArray<double_complext> & x , double y);
template void operator-=<double> (cuNDArray<double_complext> & x , double y);
template void operator*=<double> (cuNDArray<double_complext> & x , double y);
template void operator/=<double> (cuNDArray<double_complext> & x , double y);


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

template void operator+=<float_complext> (cuNDArray<float_complext> & x , cuNDArray<float_complext> & y);
template void operator-=<float_complext> (cuNDArray<float_complext> & x , cuNDArray<float_complext> & y);
template void operator*=<float_complext> (cuNDArray<float_complext> & x , cuNDArray<float_complext> & y);
template void operator/=<float_complext> (cuNDArray<float_complext> & x , cuNDArray<float_complext> & y);

template void operator+=<float_complext> (cuNDArray<float_complext> & x , float_complext y);
template void operator-=<float_complext> (cuNDArray<float_complext> & x , float_complext y);
template void operator*=<float_complext> (cuNDArray<float_complext> & x , float_complext y);
template void operator/=<float_complext> (cuNDArray<float_complext> & x , float_complext y);



template EXPORTGPUCORE class cuNDArray<double_complext>;

template void operator+=<double_complext> (cuNDArray<double_complext> & x , cuNDArray<double_complext> & y);
template void operator-=<double_complext> (cuNDArray<double_complext> & x , cuNDArray<double_complext> & y);
template void operator*=<double_complext> (cuNDArray<double_complext> & x , cuNDArray<double_complext> & y);
template void operator/=<double_complext> (cuNDArray<double_complext> & x , cuNDArray<double_complext> & y);

template void operator+=<double_complext> (cuNDArray<double_complext> & x , double_complext y);
template void operator-=<double_complext> (cuNDArray<double_complext> & x , double_complext y);
template void operator*=<double_complext> (cuNDArray<double_complext> & x , double_complext y);
template void operator/=<double_complext> (cuNDArray<double_complext> & x , double_complext y);



