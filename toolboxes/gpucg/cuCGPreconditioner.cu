#include "cuCGPreconditioner.h"

//#include <cuComplex.h>
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

template <class T> int cuCGPrecondWeight<T>::apply(cuNDArray<T>* in, cuNDArray<T>* out)
{
  if (in->get_number_of_elements() != out->get_number_of_elements()) {
    std::cerr << "cuCGPreconWeight::apply : input and output dimensions mismatch" << std::endl;
    return -1;
  }

  if (in->get_number_of_elements() % weights_->get_number_of_elements()) {
    std::cerr << "cuCGPreconWeight::apply : input dimensions don't match weights dimensions" << std::endl;
    return -1;
  }

  unsigned int num_frames = in->get_number_of_elements() / weights_->get_number_of_elements();

  dim3 blockDim(256);
  dim3 gridDim((unsigned int) ceil((double)weights_->get_number_of_elements()/blockDim.x), num_frames );
  weight_multiplication<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(),
						  weights_->get_data_ptr(), weights_->get_number_of_elements());

  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cuCGPreconWeight::apply : Unable to apply weights: " << 
      cudaGetErrorString(err) << std::endl;
    return -2;
  }
  
  return 0;
}


//
// Instantiation
//

template class EXPORTGPUCG cuCGPreconditioner<float>;
template class EXPORTGPUCG cuCGPreconditioner<float_complext::Type>;

template class EXPORTGPUCG cuCGPreconditioner<double>;
template class EXPORTGPUCG cuCGPreconditioner<double_complext::Type>;

template class EXPORTGPUCG cuCGPrecondWeight<float>;
template class EXPORTGPUCG cuCGPrecondWeight<float_complext::Type>;

template class EXPORTGPUCG cuCGPrecondWeight<double>;
template class EXPORTGPUCG cuCGPrecondWeight<double_complext::Type>;
