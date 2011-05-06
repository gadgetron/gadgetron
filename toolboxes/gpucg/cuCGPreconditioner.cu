#include "cuCGPreconditioner.h"

//#include <cuComplex.h>
#include "vector_td_utilities.h"

template<class T> __global__ void 
weight_multiplication( T* in, T* out, T* weight, unsigned long elements )
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < elements) {
    out[idx_in] = mul<T,T>(in[idx_in], weight[idx_in]);
  }
}

template <class T> int cuCGPrecondWeight<T>::apply(cuNDArray<T>* in, cuNDArray<T>* out)
{
  if (in->get_number_of_elements() != out->get_number_of_elements()) {
    std::cerr << "cuCGPreconWeight::apply : input and output dimensions mismatch" << std::endl;
    return -1;
  }

  if (in->get_number_of_elements() != weights_.get_number_of_elements()) {
    std::cerr << "cuCGPreconWeight::apply : input dimensions don't match weights dimensions" << std::endl;
    return -1;
  }

  dim3 blockDim(256,1,1);
  dim3 gridDim((unsigned int) ceil((double)in->get_number_of_elements()/blockDim.x), 1, 1 );
  weight_multiplication<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(),
						  weights_.get_data_ptr(), in->get_number_of_elements());

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

template class cuCGPrecondWeight<float>;
//template class cuCGPrecondWeight<float2>;
template class cuCGPrecondWeight<float_complext::Type>;
