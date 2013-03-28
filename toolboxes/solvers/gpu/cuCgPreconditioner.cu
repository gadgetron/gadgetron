#include "cuCgPreconditioner.h"
#include "vector_td_utilities.h"

#include <sstream>

namespace Gadgetron{

  template<class T> __global__ void 
  weight_multiplication( T* in, T* out, T* weight, unsigned long elements )
  {
    unsigned long idx = blockIdx.x*blockDim.x+threadIdx.x;
    if (idx < elements) {
      const unsigned int frame_offset = blockIdx.y*elements;
      out[idx+frame_offset] = in[idx+frame_offset]*weight[idx];
    }
  }

  template <class T> void
  cuCgPreconditioner<T>::set_weights( boost::shared_ptr< cuNDArray<T> > w ) 
  {
    weights_ = w;
  }
  
  template <class T> void
  cuCgPreconditioner<T>::apply(cuNDArray<T>* in, cuNDArray<T>* out)
  {
    if( !weights_.get() ){
      throw std::runtime_error( "cuCGPreconWeight::apply: weights not set");
    }
 
    if ( !in || !out || in->get_number_of_elements() != out->get_number_of_elements()) {
      throw std::runtime_error("cuCGPreconWeight::apply: input and output dimensions mismatch");
    }

    if (in->get_number_of_elements() % weights_->get_number_of_elements()) {
      throw std::runtime_error( "cuCGPreconWeight::apply: input dimensions don't match weights dimensions" );
    }
  
    cuNDArray<T> *in_int, *out_int;
    in_int = in;
    out_int = out;
  
    unsigned int num_frames = in->get_number_of_elements() / weights_->get_number_of_elements();

    dim3 blockDim(256);
    dim3 gridDim((weights_->get_number_of_elements()+blockDim.x-1)/blockDim.x, num_frames );
    weight_multiplication<<< gridDim, blockDim >>>( in_int->get_data_ptr(), out_int->get_data_ptr(),
						    weights_->get_data_ptr(), weights_->get_number_of_elements());

    cudaError_t err = cudaGetLastError();
    if( err != cudaSuccess ){
      std::stringstream ss;
      ss << "cuCGPreconWeight::apply: Unable to apply weights: " <<
	cudaGetErrorString(err);
      throw std::runtime_error(ss.str());
    }  
  }

  //
  // Instantiation
  //

  template class EXPORTGPUSOLVERS cuCgPreconditioner<float>;
  template class EXPORTGPUSOLVERS cuCgPreconditioner<float_complext>;

  template class EXPORTGPUSOLVERS cuCgPreconditioner<double>;
  template class EXPORTGPUSOLVERS cuCgPreconditioner<double_complext>;
}
