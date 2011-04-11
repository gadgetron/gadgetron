#include "cgOperatorSense.h"

__global__ void clear_array(float2* in, unsigned long int elements)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < elements) {
    in[idx_in].x = 0.0;
    in[idx_in].y = 0.0;
  }
}

__global__ void coil_sensitivity_multiplication(float2* in, float2* out, float2* csm, 
						unsigned long image_elements,
						unsigned int coils)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < image_elements) {
    for (unsigned int i = 0; i < coils; i++) {
      out[idx_in + i*image_elements] = cuCmulf(in[idx_in],csm[idx_in + i*image_elements]);
    }
  }
}

__global__ void coil_sensitivity_conj_mult_sum(float2* in, float2* out, float2* csm, 
					       unsigned long image_elements,
					       unsigned int coils)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < image_elements) {
    for (unsigned int i = 0; i < coils; i++) {
      float2 tmp = cuCmulf(in[idx_in + i*image_elements],cuConjf(csm[idx_in + i*image_elements]));
      out[idx_in].x += tmp.x;
      out[idx_in].y += tmp.y; 
    }
  }
}


int cgOperatorSense::mult_csm(cuNDArray<float2>* in, cuNDArray<float2>* out)
{

  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)in->get_number_of_elements()/blockDim.x), 1, 1 );
  coil_sensitivity_multiplication<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(),
							    csm_->get_data_ptr(), in->get_number_of_elements(), coils_);

  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cgOperatorSense::mult_csm Unable to multiply with coil sensitivities: " << 
      cudaGetErrorString(err) << std::endl;
    return -1;
  }
  return 0;
}


int cgOperatorSense:: mult_csm_conj_sum(cuNDArray<float2>* in, cuNDArray<float2>* out)
{

  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)out->get_number_of_elements()/blockDim.x), 1, 1 );
  coil_sensitivity_conj_mult_sum<<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(),
							   csm_->get_data_ptr(),out->get_number_of_elements(),
							   coils_);
  
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cgOperatorCartesianSense::mult_EM : Unable to combine coils " << 
      cudaGetErrorString(err) << std::endl;
    return -1;
  }

  return 0;
}


int cgOperatorSense::clear(cuNDArray<float2>* in)
{
  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)in->get_number_of_elements()/blockDim.x), 1, 1 );

  clear_array<<< gridDim, blockDim >>>( in->get_data_ptr(), in->get_number_of_elements());

  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cgOperatorCartesianSense::clear : Error during kernel call: " << cudaGetErrorString(err) << std::endl;
    return -1;
  }

  return 0;
}


int cgOperatorSense::mult_MH_M(cuNDArray<float2>* in, cuNDArray<float2>* out, bool accumulate)
{
  cuNDArray<float2> tmp;
  if (!tmp.create(dimensions_out_)) {
    std::cerr << "cgOperatorCartesianSense::mult_MH_M: Unable to create temporary storage" << std::endl;
    return -1;
  }

  if (mult_M(in, &tmp, false) < 0) {
    std::cerr << "cgOperatorCartesianSense::mult_MH_M: Unable to perform mult_M" << std::endl;
    return -2;
  }

  if (mult_MH(&tmp, out, accumulate) < 0) {
    std::cerr << "cgOperatorCartesianSense::mult_MH_M: Unable to perform mult_M" << std::endl;
    return -2;
  }

  return 0;
}
