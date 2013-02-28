#include "cuNDArray.h"
#include "vector_td.h"
#include <thrust/iterator/constant_iterator.h>
#include <stdexcept>
#include <sstream>
#include "vector_td_utilities.h"
#include "check_CUDA.h"


using namespace Gadgetron;

template<typename T>
class cuNDA_modulus : public thrust::unary_function<T,T>
{
public:
	cuNDA_modulus(int x):mod(x) {};
	 __host__ __device__ T operator()(const T &y) const {return y%mod;}
private:
	const int mod;
};

template<class T,class S> static bool compatible_dimensions(cuNDArray<T> & x, cuNDArray<S> & y){
	bool retVal = true;
	for (int i = 0; i < y.get_number_of_dimensions(); i++){
		retVal &= (x.get_size(i) == y.get_size(i));
	}
	return retVal;
}
template<class T,class S, class F>  void equals_transform(cuNDArray<T> & x,cuNDArray<S> & y){
	if (x.dimensions_equal(&y)){
			thrust::transform(x.begin(), x.end(), y.begin(), x.begin(), F());
		} else if (compatible_dimensions(x,y))
		{
			typedef thrust::transform_iterator<cuNDA_modulus<int>,thrust::counting_iterator<int>, int> transform_it;
			transform_it indices = thrust::make_transform_iterator(thrust::make_counting_iterator(0),cuNDA_modulus<int>(y.get_number_of_elements()));
			thrust::permutation_iterator<thrust::device_ptr<S>,transform_it> p = thrust::make_permutation_iterator(y.begin(),indices);
			thrust::transform(x.begin(),x.end(),p,x.begin(),F());
		} else {
			BOOST_THROW_EXCEPTION(runtime_error("cuNDArrays have incompatible dimensions"));
		}
}

template <class T> 
Gadgetron::cuNDArray<T>::cuNDArray() : NDArray<T>::NDArray()
{ 
  cudaGetDevice(&this->device_); 
}

template <class T> 
Gadgetron::cuNDArray<T>::cuNDArray(std::vector<unsigned int> *dimensions) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions);
}

template <class T>
Gadgetron::cuNDArray<T>::cuNDArray(std::vector<unsigned int> *dimensions, int device_no) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions,device_no);
}

template <class T>
Gadgetron::cuNDArray<T>::cuNDArray(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions,data,delete_data_on_destruct);
}


template <class T>
Gadgetron::cuNDArray<T>::cuNDArray(boost::shared_ptr<std::vector<unsigned int> > dimensions) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions.get());
}

template <class T>
Gadgetron::cuNDArray<T>::cuNDArray(boost::shared_ptr<std::vector<unsigned int> > dimensions, int device_no) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions.get(),device_no);
}

template <class T>
Gadgetron::cuNDArray<T>::cuNDArray(boost::shared_ptr<std::vector<unsigned int> > dimensions, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  create(dimensions.get(),data,delete_data_on_destruct);
}

template <class T>
Gadgetron::cuNDArray<T>::cuNDArray(hoNDArray<T> *a) : NDArray<T>::NDArray()
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
Gadgetron::cuNDArray<T>::cuNDArray(const hoNDArray<T>& a) : NDArray<T>::NDArray()
{
  cudaGetDevice(&this->device_);
  this->dimensions_ = a.get_dimensions();

  allocate_memory();
	if (cudaMemcpy(this->data_, a.get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
	  deallocate_memory();
	  this->data_ = 0;
	  this->dimensions_->clear();

  }
}

template <class T> 
Gadgetron::cuNDArray<T>::cuNDArray(const cuNDArray<T>& a)
{
  cudaGetDevice(&this->device_);
  this->data_ = 0;
  this->dimensions_ = a.dimensions_;

  allocate_memory();
	if (a.device_ == this->device_) {
	  CUDA_CALL(cudaMemcpy(this->data_, a.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice));
	} else {
	  //This memory is on a different device, we must move it.
	  cudaSetDevice(a.device_);
	  boost::shared_ptr< hoNDArray<T> > tmp = a.to_host();
	  cudaSetDevice(this->device_);

	  cudaError_t err = cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice);
	  if (err !=cudaSuccess) {
			deallocate_memory();
			this->data_ = 0;
			this->dimensions_->clear();
			BOOST_THROW_EXCEPTION(cuda_error(err));
	  }
	}

}

template <class T> 
Gadgetron::cuNDArray<T>& Gadgetron::cuNDArray<T>::operator=(const cuNDArray<T>& rhs)
{
  int cur_device; 
  CUDA_CALL(cudaGetDevice(&cur_device));
  
  bool dimensions_match = this->dimensions_equal(&rhs);
  
  if (dimensions_match && (rhs.device_ == cur_device) && (cur_device == this->device_)) {
    
    CUDA_CALL(cudaMemcpy(this->data_, rhs.data_, this->elements_*sizeof(T), cudaMemcpyDeviceToDevice));
  } 
  else {
    
    CUDA_CALL(cudaSetDevice(this->device_));
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
	BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::operator=: failed to copy data (2)"));

      }
    } else {

      if( cudaSetDevice(rhs.device_) != cudaSuccess) {

	cudaSetDevice(cur_device);
	BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::operator=: unable to set device no (2)"));

      }

      boost::shared_ptr< hoNDArray<T> > tmp = rhs.to_host();
      if( cudaSetDevice(this->device_) != cudaSuccess) {

	cudaSetDevice(cur_device);
	BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::operator=: unable to set device no (3)"));

      }
      
      if (cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {

	cudaSetDevice(cur_device);
	BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::operator=: failed to copy data (3)"));

      }
    }

    if( cudaSetDevice(cur_device) != cudaSuccess) {
      BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::operator=: unable to restore to current device"));

    }
  }
  
  return *this;
}

template <class T>
Gadgetron::cuNDArray<T>::~cuNDArray()
{ 
  if (this->delete_data_on_destruct_) 
    deallocate_memory();  
}

template <class T>
void Gadgetron::cuNDArray<T>::create(std::vector<unsigned int> *dimensions)
{
  NDArray<T>::create(dimensions);
}

template <class T>
void Gadgetron::cuNDArray<T>::create(std::vector<unsigned int> *dimensions, int device_no)
{
  if (device_no < 0){
    BOOST_THROW_EXCEPTION(cuda_error("cuNDArray::create: illegal device no"));

  }
  
  this->device_ = device_no; 
  NDArray<T>::create(dimensions);
}

template <class T>
void Gadgetron::cuNDArray<T>::create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct)
{
  if (!data) {
    BOOST_THROW_EXCEPTION(runtime_error("cuNDArray::create: 0x0 pointer provided"));
  }
  
  int tmp_device; 
  if( cudaGetDevice(&tmp_device) != cudaSuccess) {
    BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::create: Unable to query for device"));

  }
  
  cudaDeviceProp deviceProp; 
  if( cudaGetDeviceProperties( &deviceProp, tmp_device) != cudaSuccess) {
    BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::create: Unable to query device properties"));

  }
  
  if (deviceProp.unifiedAddressing) {
    cudaPointerAttributes attrib;
    if (cudaPointerGetAttributes(&attrib, data) != cudaSuccess) {
      BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::create: Unable to determine attributes of pointer"));

    }
    this->device_ = attrib.device;
  } else {
    this->device_ = tmp_device;
  }
  
  NDArray<T>::create(dimensions, data, delete_data_on_destruct);
}

template <class T>
boost::shared_ptr< Gadgetron::cuNDArray<T> > Gadgetron::cuNDArray<T>::allocate(std::vector<unsigned int> *dimensions)
{
  int tmp_device; cudaGetDevice(&tmp_device);
  return allocate( dimensions, tmp_device );
}

template <class T>
boost::shared_ptr< Gadgetron::cuNDArray<T> > Gadgetron::cuNDArray<T>::allocate(std::vector<unsigned int> *dimensions, int device_no)
{
  boost::shared_ptr< cuNDArray<T> > ret( new cuNDArray<T> );
  
  ret->create(dimensions, device_no);
  
  return ret;
}

template <class T>
boost::shared_ptr< Gadgetron::hoNDArray<T> > Gadgetron::cuNDArray<T>::to_host() const
{
  boost::shared_ptr< hoNDArray<T> > ret = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>);
  ret->create(this->dimensions_.get());
    if (cudaMemcpy(ret->get_data_ptr(), this->data_, this->elements_*sizeof(T), cudaMemcpyDeviceToHost) != cudaSuccess) {
      BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::to_host(): failed to copy memory from device"));
    }

  return ret;
}

template <class T> 
void Gadgetron::cuNDArray<T>::permute(std::vector<unsigned int> *dim_order, NDArray<T> *out, int shift_mode)
{

  cuNDArray<T>* out_int = 0x0;
  
  //Check ordering array
  if (dim_order->size() > this->dimensions_->size()) {
    BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::permute - Invalid length of dimension ordering array"));
  }
  
  std::vector<unsigned int> dim_count(this->dimensions_->size(),0);
  for (unsigned int i = 0; i < dim_order->size(); i++) {
    if ((*dim_order)[i] >= this->dimensions_->size()) {
      BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::permute - Invalid dimension order array"));
    }
    dim_count[(*dim_order)[i]]++;
  }

  //Create an internal array to store the dimensions
  std::vector<unsigned int> dim_order_int;
  
  //Check that there are no duplicate dimensions
  for (unsigned int i = 0; i < dim_order->size(); i++) {
    if (dim_count[(*dim_order)[i]] != 1) {
    	BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::permute - Invalid dimension order array (duplicates)"));
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
      BOOST_THROW_EXCEPTION(runtime_error("cuNDArray::permute: failed to dynamic cast out array pointer"));
    }
    for (unsigned int i = 0; i < dim_order_int.size(); i++) {
      if ((*this->dimensions_)[dim_order_int[i]] != out_int->get_size(i)) {
	BOOST_THROW_EXCEPTION(runtime_error("cuNDArray::permute: Dimensions of output array do not match the input array"));
      }
    }
  }
  
  cuNDArray_permute(this, out_int, &dim_order_int, shift_mode);
}

template <class T> 
void Gadgetron::cuNDArray<T>::set_device(int device)
{
  if( device_ == device )
    return;

  int cur_device; 
  if( cudaGetDevice(&cur_device) != cudaSuccess) {
    BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::set_device: unable to get device no"));

  }

  if( cur_device != device_ && cudaSetDevice(device_) != cudaSuccess) {
    BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::set_device: unable to set device no"));

  }
  
  boost::shared_ptr< hoNDArray<T> > tmp = to_host();

  deallocate_memory();

  
  if( cudaSetDevice(device) != cudaSuccess) {
	  cudaSetDevice(cur_device);
	  BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::set_device: unable to set device no (2)"));
  }
 
  device_ = device;
  
  allocate_memory();
  if (cudaMemcpy(this->data_, tmp->get_data_ptr(), this->elements_*sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
	  cudaSetDevice(cur_device);
    BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::set_device: failed to copy data"));
  }

  if( cudaSetDevice(cur_device) != cudaSuccess) {
    BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::set_device: unable to restore device to current device"));
  }

}

template <class T>
void Gadgetron::cuNDArray<T>::allocate_memory()
{
  deallocate_memory();
  
  this->elements_ = 1;
  if (this->dimensions_->empty())
  	BOOST_THROW_EXCEPTION( runtime_error("cuNDArray::allocate_memory() : dimensions is empty."));
  for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
    this->elements_ *= (*this->dimensions_)[i];
  } 
  
  size_t size = this->elements_ * sizeof(T);
  
  int device_no_old;
  if (cudaGetDevice(&device_no_old) != cudaSuccess) {
    BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::allocate_memory: unable to get device no"));
  }
  
  if (device_ != device_no_old) {
    if (cudaSetDevice(device_) != cudaSuccess) {
      BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::allocate_memory: unable to set device no"));
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
    BOOST_THROW_EXCEPTION(bad_alloc(err.str()));

  }
  
  if (device_ != device_no_old) {
    if (cudaSetDevice(device_no_old) != cudaSuccess) {
      BOOST_THROW_EXCEPTION( cuda_error("cuNDArray::allocate_memory: unable to restore device no"));
    }
  }

}

template <class T>
void Gadgetron::cuNDArray<T>::deallocate_memory()
{
  if (this->data_) {
    
    int device_no_old;
    CUDA_CALL(cudaGetDevice(&device_no_old));
    if (device_ != device_no_old) {
      CUDA_CALL(cudaSetDevice(device_));
    }
    
    CUDA_CALL(cudaFree(this->data_));
    if (device_ != device_no_old) {
      CUDA_CALL(cudaSetDevice(device_no_old));
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
struct cuNDA_abs<Gadgetron::vector_td<T,D> > : public thrust::unary_function<Gadgetron::vector_td<T,D>, Gadgetron::vector_td<T,D> >
{
 __host__ __device__ Gadgetron::vector_td<T,D> operator()(const Gadgetron::vector_td<T,D> &x) const {
	 Gadgetron::vector_td<T,D> res;
	 for (int i = 0; i < D; i++) res[i] = abs(x[i]);
	 return res;
 	 }
};

template <class T>
void Gadgetron::cuNDArray<T>::abs()
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
struct cuNDA_sqrt<Gadgetron::vector_td<T,D> > : public thrust::unary_function<Gadgetron::vector_td<T,D>, Gadgetron::vector_td<T,D> >
{
 __host__ __device__ Gadgetron::vector_td<T,D> operator()(const Gadgetron::vector_td<T,D> &x) const {
	 Gadgetron::vector_td<T,D> res;
	 for (int i = 0; i < D; i++) res[i] = sqrt(x[i]);
	 return res;
 	 }
};

template <class T>
void Gadgetron::cuNDArray<T>::sqrt()
{
	thrust::device_ptr<T> devPtr = this->get_device_ptr();
	thrust::transform(devPtr,devPtr+this->get_number_of_elements(),devPtr,cuNDA_sqrt<T>());
}



template <class T>
void Gadgetron::cuNDArray<T>::clear()
{
	cudaMemset(this->get_data_ptr(),0,sizeof(T)*this->get_number_of_elements());
}

template <class T>
void Gadgetron::cuNDArray<T>::fill(T val)
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
struct cuNDA_reciprocal<Gadgetron::vector_td<T,D> > : public thrust::unary_function<Gadgetron::vector_td<T,D>, Gadgetron::vector_td<T,D> >
{
 __host__ __device__ Gadgetron::vector_td<T,D> operator()(const Gadgetron::vector_td<T,D> &x) const {
	 Gadgetron::vector_td<T,D> res;
	 for (int i = 0; i < D; i++) res[i] = 1/(x[i]);
	 return res;
 	 }
};
template <class T>
void Gadgetron::cuNDArray<T>::reciprocal()
{
	thrust::device_ptr<T> devPtr = this->get_device_ptr();
	thrust::transform(devPtr,devPtr+this->get_number_of_elements(),devPtr,cuNDA_reciprocal<T>());
}

template<typename T>
struct cuNDA_reciprocal_sqrt : public thrust::unary_function<T,T>
{
 __host__ __device__ T operator()(const T &x) const {return 1/sqrt(x);}
};

template<typename T, unsigned int D>
struct cuNDA_reciprocal_sqrt<Gadgetron::vector_td<T,D> > : public thrust::unary_function<Gadgetron::vector_td<T,D>, Gadgetron::vector_td<T,D> >
{
 __host__ __device__ Gadgetron::vector_td<T,D> operator()(const Gadgetron::vector_td<T,D> &x) const {
	 Gadgetron::vector_td<T,D> res;
	 for (int i = 0; i < D; i++) res[i] = 1/sqrt((x[i]));
	 return res;
 	 }
};
template <class T>
void Gadgetron::cuNDArray<T>::reciprocal_sqrt()
{
	thrust::device_ptr<T> devPtr = this->get_device_ptr();
	thrust::transform(devPtr,devPtr+this->get_number_of_elements(),devPtr,cuNDA_reciprocal_sqrt<T>());
}


template<class T> void Gadgetron::operator+= (Gadgetron::cuNDArray<T> & x , Gadgetron::cuNDArray<T> & y){
	equals_transform< T,T,thrust::plus<T> >(x,y);
}

template<class T> void Gadgetron::operator+= (Gadgetron::cuNDArray<T> & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::plus<T>());
}

template<class T> void Gadgetron::operator*= (Gadgetron::cuNDArray<T> & x , Gadgetron::cuNDArray<T> & y){
	equals_transform< T,T,thrust::multiplies<T> >(x,y);
}

template<class T> void Gadgetron::operator*= (Gadgetron::cuNDArray<T> & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::multiplies<T>());
}

template<class T> void Gadgetron::operator-= (Gadgetron::cuNDArray<T> & x , Gadgetron::cuNDArray<T> & y){
	equals_transform< T,T,thrust::minus<T> >(x,y);
}

template<class T> void Gadgetron::operator-= (Gadgetron::cuNDArray<T> & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), thrust::minus<T>());
}
template<class T> void Gadgetron::operator/= (Gadgetron::cuNDArray<T> & x , Gadgetron::cuNDArray<T> & y){
	equals_transform< T,T,thrust::divides<T> >(x,y);
}

template<class T> void Gadgetron::operator/= (Gadgetron::cuNDArray<T> & x , T y){
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




template<class T> void Gadgetron::operator+= (Gadgetron::cuNDArray< complext<T> > & x , Gadgetron::cuNDArray<T> & y){
	equals_transform< complext<T>,T,cuNDA_plus<T> >(x,y);

}

template<class T> void Gadgetron::operator+= (Gadgetron::cuNDArray<complext<T> > & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_plus<T>());
}

template<class T> void Gadgetron::operator*= (Gadgetron::cuNDArray< complext<T> > & x , Gadgetron::cuNDArray<T> & y){
	equals_transform< complext<T>,T,cuNDA_multiplies<T> >(x,y);
}

template<class T> void Gadgetron::operator*= (Gadgetron::cuNDArray<complext<T> > & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_multiplies<T>());
}

template<class T> void Gadgetron::operator-= (Gadgetron::cuNDArray< complext<T> > & x , Gadgetron::cuNDArray<T> & y){
	equals_transform< complext<T>,T,cuNDA_minus<T> >(x,y);
}

template<class T> void Gadgetron::operator-= (Gadgetron::cuNDArray<complext<T> > & x , T y){
	thrust::constant_iterator<T> iter(y);
	thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_minus<T>());
}
template<class T> void Gadgetron::operator/= (Gadgetron::cuNDArray< complext<T> > & x , Gadgetron::cuNDArray<T> & y){
	equals_transform< complext<T>,T,cuNDA_divides<T> >(x,y);
}

template<class T> void Gadgetron::operator/= (Gadgetron::cuNDArray<complext<T> > & x , T y){
	thrust::constant_iterator<T> iter(y);

	thrust::transform(x.begin(), x.end(), iter, x.begin(), cuNDA_divides<T>());
}


//
// Instantiation
//

template EXPORTGPUCORE class Gadgetron::cuNDArray< int >;
/*
template EXPORTGPUCORE class cuNDArray< int2 >;
template EXPORTGPUCORE class cuNDArray< int3 >;
template EXPORTGPUCORE class cuNDArray< int4 >;
*/
template void Gadgetron::operator+=<int> (cuNDArray<int> & x , cuNDArray<int> & y);
template void Gadgetron::operator-=<int> (cuNDArray<int> & x , cuNDArray<int> & y);
template void Gadgetron::operator*=<int> (cuNDArray<int> & x , cuNDArray<int> & y);
template void Gadgetron::operator/=<int> (cuNDArray<int> & x , cuNDArray<int> & y);

template void Gadgetron::operator+=<int> (cuNDArray<int> & x , int y);
template void Gadgetron::operator-=<int> (cuNDArray<int> & x , int y);
template void Gadgetron::operator*=<int> (cuNDArray<int> & x , int y);
template void Gadgetron::operator/=<int> (cuNDArray<int> & x , int y);

template EXPORTGPUCORE class Gadgetron::cuNDArray< unsigned int >;

/*
template EXPORTGPUCORE class Gadgetron::cuNDArray< uint2 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< uint3 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< uint4 >;
*/

template void Gadgetron::operator+=<unsigned int> (cuNDArray<unsigned int> & x , cuNDArray<unsigned int> & y);
template void Gadgetron::operator-=<unsigned int> (cuNDArray<unsigned int> & x , cuNDArray<unsigned int> & y);
template void Gadgetron::operator*=<unsigned int> (cuNDArray<unsigned int> & x , cuNDArray<unsigned int> & y);
template void Gadgetron::operator/=<unsigned int> (cuNDArray<unsigned int> & x , cuNDArray<unsigned int> & y);

template void Gadgetron::operator+=<unsigned int> (cuNDArray<unsigned int> & x , unsigned int y);
template void Gadgetron::operator-=<unsigned int> (cuNDArray<unsigned int> & x , unsigned int y);
template void Gadgetron::operator*=<unsigned int> (cuNDArray<unsigned int> & x , unsigned int y);
template void Gadgetron::operator/=<unsigned int> (cuNDArray<unsigned int> & x , unsigned int y);

template EXPORTGPUCORE class Gadgetron::cuNDArray< float >;



/*template EXPORTGPUCORE class Gadgetron::cuNDArray< float2 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< float3 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< float4 >;
*/

template void Gadgetron::operator+=<float> (Gadgetron::cuNDArray<float> & x , Gadgetron::cuNDArray<float> & y);
template void Gadgetron::operator-=<float> (Gadgetron::cuNDArray<float> & x , Gadgetron::cuNDArray<float> & y);
template void Gadgetron::operator*=<float> (Gadgetron::cuNDArray<float> & x , Gadgetron::cuNDArray<float> & y);
template void Gadgetron::operator/=<float> (Gadgetron::cuNDArray<float> & x , Gadgetron::cuNDArray<float> & y);

template void Gadgetron::operator+=<float> (Gadgetron::cuNDArray<float> & x , float y);
template void Gadgetron::operator-=<float> (Gadgetron::cuNDArray<float> & x , float y);
template void Gadgetron::operator*=<float> (Gadgetron::cuNDArray<float> & x , float y);
template void Gadgetron::operator/=<float> (Gadgetron::cuNDArray<float> & x , float y);


template void Gadgetron::operator+=<float> (cuNDArray<float_complext> & x , cuNDArray<float> & y);
template void Gadgetron::operator-=<float> (cuNDArray<float_complext> & x , cuNDArray<float> & y);
template void Gadgetron::operator*=<float> (cuNDArray<float_complext> & x , cuNDArray<float> & y);
template void Gadgetron::operator/=<float> (cuNDArray<float_complext> & x , cuNDArray<float> & y);


template void Gadgetron::operator+=<float> (cuNDArray<float_complext> & x , float y);
template void Gadgetron::operator-=<float> (cuNDArray<float_complext> & x , float y);
template void Gadgetron::operator*=<float> (cuNDArray<float_complext> & x , float y);
template void Gadgetron::operator/=<float> (cuNDArray<float_complext> & x , float y);


template EXPORTGPUCORE class Gadgetron::cuNDArray< double >;

/*
template EXPORTGPUCORE class Gadgetron::cuNDArray< double2 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< double3 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< double4 >;
*/

template void Gadgetron::operator+=<double> (Gadgetron::cuNDArray<double> & x , Gadgetron::cuNDArray<double> & y);
template void Gadgetron::operator-=<double> (Gadgetron::cuNDArray<double> & x , Gadgetron::cuNDArray<double> & y);
template void Gadgetron::operator*=<double> (Gadgetron::cuNDArray<double> & x , Gadgetron::cuNDArray<double> & y);
template void Gadgetron::operator/=<double> (Gadgetron::cuNDArray<double> & x , Gadgetron::cuNDArray<double> & y);

template void Gadgetron::operator+=<double> (Gadgetron::cuNDArray<double> & x , double y);
template void Gadgetron::operator-=<double> (Gadgetron::cuNDArray<double> & x , double y);
template void Gadgetron::operator*=<double> (Gadgetron::cuNDArray<double> & x , double y);
template void Gadgetron::operator/=<double> (Gadgetron::cuNDArray<double> & x , double y);


template EXPORTGPUCORE void Gadgetron::operator+=<double> (cuNDArray<double_complext> & x , cuNDArray<double> & y);
template EXPORTGPUCORE void Gadgetron::operator-=<double> (cuNDArray<double_complext> & x , cuNDArray<double> & y);
template EXPORTGPUCORE void Gadgetron::operator*=<double> (cuNDArray<double_complext> & x , cuNDArray<double> & y);
template EXPORTGPUCORE void Gadgetron::operator /=<double> (cuNDArray<double_complext> & x , cuNDArray<double> & y);

template void Gadgetron::operator+=<double> (cuNDArray<double_complext> & x , double y);
template void Gadgetron::operator-=<double> (cuNDArray<double_complext> & x , double y);
template void Gadgetron::operator*=<double> (cuNDArray<double_complext> & x , double y);
template void Gadgetron::operator/=<double> (cuNDArray<double_complext> & x , double y);


template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::intd1 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::intd2 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::intd3 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::intd4 >;

template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::uintd1 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::uintd2 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::uintd3 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::uintd4 >;

template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::floatd1 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::floatd2 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::floatd3 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::floatd4 >;

template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::doubled1 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::doubled2 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::doubled3 >;
template EXPORTGPUCORE class Gadgetron::cuNDArray< Gadgetron::doubled4 >;

template EXPORTGPUCORE class Gadgetron::cuNDArray<float_complext>;

template void Gadgetron::operator+=<float_complext> (cuNDArray<float_complext> & x , cuNDArray<float_complext> & y);
template void Gadgetron::operator-=<float_complext> (cuNDArray<float_complext> & x , cuNDArray<float_complext> & y);
template void Gadgetron::operator*=<float_complext> (cuNDArray<float_complext> & x , cuNDArray<float_complext> & y);
template void Gadgetron::operator/=<float_complext> (cuNDArray<float_complext> & x , cuNDArray<float_complext> & y);

template void Gadgetron::operator+=<float_complext> (cuNDArray<float_complext> & x , float_complext y);
template void Gadgetron::operator-=<float_complext> (cuNDArray<float_complext> & x , float_complext y);
template void Gadgetron::operator*=<float_complext> (cuNDArray<float_complext> & x , float_complext y);
template void Gadgetron::operator/=<float_complext> (cuNDArray<float_complext> & x , float_complext y);



template EXPORTGPUCORE class Gadgetron::cuNDArray<double_complext>;

template void Gadgetron::operator+=<double_complext> (cuNDArray<double_complext> & x , cuNDArray<double_complext> & y);
template void Gadgetron::operator-=<double_complext> (cuNDArray<double_complext> & x , cuNDArray<double_complext> & y);
template void Gadgetron::operator*=<double_complext> (cuNDArray<double_complext> & x , cuNDArray<double_complext> & y);
template void Gadgetron::operator/=<double_complext> (cuNDArray<double_complext> & x , cuNDArray<double_complext> & y);

template void Gadgetron::operator+=<double_complext> (cuNDArray<double_complext> & x , double_complext y);
template void Gadgetron::operator-=<double_complext> (cuNDArray<double_complext> & x , double_complext y);
template void Gadgetron::operator*=<double_complext> (cuNDArray<double_complext> & x , double_complext y);
template void Gadgetron::operator/=<double_complext> (cuNDArray<double_complext> & x , double_complext y);


