#include "cuCGPrecondWeights.h"
#include "vector_td_utilities.h"

template<class T> __global__ void 
weight_multiplication( T* in, T* out, T* weight, unsigned long elements )
{
  unsigned long idx = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx < elements) {
    const unsigned int frame_offset = blockIdx.y*elements;
    out[idx+frame_offset] = mul<T>(in[idx+frame_offset], weight[idx]);
  }
}

template <class T> int 
cuCGPrecondWeights<T>::set_weights( boost::shared_ptr< cuNDArray<T> > w ) 
{
  int ret1 = 0, ret2 = 0;
  if( w->get_device() != this->device_ ){
    ret1 = this->set_device();
    if( ret1 == 0 ) weights_ = boost::shared_ptr< cuNDArray<T> >(new cuNDArray<T>(*w.get()));
    ret2 = this->restore_device();
  }
  else    
    weights_ = w;
  
  if( ret1 == 0 && ret2 == 0 )
    return 0;
  else{
    std::cout << std::endl << "cuCGPrecondWeight::set_weights failed" << std::endl;
    return -1;
  }
}

template <class T> int 
cuCGPrecondWeights<T>::apply(cuNDArray<T>* in, cuNDArray<T>* out)
{
	if( !weights_.get() ){
    std::cerr << "cuCGPreconWeight::apply: weights not set" << std::endl;
    return -1;
  }
 
  if ( !in || in->get_number_of_elements() != out->get_number_of_elements()) {
    std::cerr << "cuCGPreconWeight::apply: input and output dimensions mismatch" << std::endl;
    return -1;
  }

  if (in->get_number_of_elements() % weights_->get_number_of_elements()) {
    std::cerr << "cuCGPreconWeight::apply: input dimensions don't match weights dimensions" << std::endl;
    return -1;
  }
  
  int ret = this->set_device();
  if( ret < 0 ){
    std::cerr << "cuCGPreconWeight::apply: unable to set device" << std::endl;
    return -1;
  }
  
  cuNDArray<T> *in_int, *out_int;

  if( this->device_ != in->get_device() )
    in_int = new cuNDArray<T>(*in);
  else 
    in_int = in;

  if( this->device_ != out->get_device() )
    out_int = new cuNDArray<T>(*out);
  else 
    out_int = out;
  
  unsigned int num_frames = in->get_number_of_elements() / weights_->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((weights_->get_number_of_elements()+blockDim.x-1)/blockDim.x, num_frames );
  weight_multiplication<<< gridDim, blockDim >>>( in_int->get_data_ptr(), out_int->get_data_ptr(),
						  weights_->get_data_ptr(), weights_->get_number_of_elements());

  if( this->device_ != in->get_device() )
    delete in_int;

  if( this->device_ != out->get_device() ){
    *out = *out_int;  
    delete out_int;
  }

  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cuCGPreconWeight::apply: Unable to apply weights: " << 
      cudaGetErrorString(err) << std::endl;
    return -2;
  }
  
  ret = this->restore_device();
  if( ret < 0 ){
    std::cerr << "cuCGPreconWeight::apply: unable to restore device" << std::endl;
    return -3;
  }
  
  return 0;
}


//
// Instantiation
//

template class EXPORTSOLVERS cuCGPrecondWeights<float>;
template class EXPORTSOLVERS cuCGPrecondWeights<float_complext::Type>;

template class EXPORTSOLVERS cuCGPrecondWeights<double>;
template class EXPORTSOLVERS cuCGPrecondWeights<double_complext::Type>;
