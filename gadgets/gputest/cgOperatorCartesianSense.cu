#include "cgOperatorCartesianSense.h"

#include "cuNDFFT.h"

__global__ void clear_array(float2* in, unsigned long int elements)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < elements) {
    in[idx_in].x = 0.0;
    in[idx_in].y = 0.0;
  }
}


__global__ void sample_array_kernel(float2* in, float2* out, 
				    unsigned int* idx, 
				    unsigned long image_elements,
				    unsigned long int samples,
				    unsigned int coils)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < samples) {
    for (unsigned int i = 0; i < coils; i++) {
      out[idx_in + i*samples].x += in[idx[idx_in] + i*image_elements].x;
      out[idx_in + i*samples].y += in[idx[idx_in] + i*image_elements].y;
    }
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

int cgOperatorCartesianSense::mult_M(cuNDArray<float2>* in, 
				     cuNDArray<float2>* out, 
				     bool accumulate)
{
  if (!(in->dimensions_equal(dimensions_)) ||
      !(out->dimensions_equal(dimensions_out_)) ) {

    std::cerr << "cgOperatorCartesianSense::mult_M dimensions mismatch" << std::endl;

    /*
    std::cout << "dimensions_ = ";
    for (unsigned int i = 0; i < dimensions_.size(); i++) {
      std::cout << " " << dimensions_[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "dimensions_out_ = ";
    for (unsigned int i = 0; i < dimensions_out_.size(); i++) {
      std::cout << " " << dimensions_out_[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "in->dimensions_ = ";
    for (unsigned int i = 0; i < in->get_dimensions().size(); i++) {
      std::cout << " " << in->get_dimensions()[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "out->dimensions_ = ";
    for (unsigned int i = 0; i < out->get_dimensions().size(); i++) {
      std::cout << " " << out->get_dimensions()[i] << " ";
    }
    std::cout << std::endl;
    */

    return -1;
  }

  cuNDArray<float2> tmp;
  std::vector<unsigned int> full_dimensions = dimensions_;
  full_dimensions.push_back(coils_);

  if (!tmp.create(full_dimensions)) {
    std::cerr << "cgOperatorCartesianSense::mult_M unable to allocate temp array" << std::endl;
    return -1;    
  }

  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)in->get_number_of_elements()/blockDim.x), 1, 1 );
  coil_sensitivity_multiplication<<< gridDim, blockDim >>>( in->get_data_ptr(), tmp.get_data_ptr(),
							    csm_->get_data_ptr(), in->get_number_of_elements(), coils_);


  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cgOperatorCartesianSense::mult_M : Unable to multiply with coil sensitivities: " << 
      cudaGetErrorString(err) << std::endl;
    return -1;
  }

  cuNDFFT ft;
  std::vector<unsigned int> ft_dims;
  for (unsigned int i = 0; i < dimensions_.size(); i++) {
    ft_dims.push_back(i);
  }

  ft.fft(&tmp, ft_dims);

  if (!accumulate) clear(out);


  gridDim = dim3((unsigned int) ceil((double)idx_->get_number_of_elements()/blockDim.x), 1, 1 );
  sample_array_kernel<<< gridDim, blockDim >>>( tmp.get_data_ptr(), out->get_data_ptr(), idx_->get_data_ptr(),
							    in->get_number_of_elements(), idx_->get_number_of_elements(), coils_);
  err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cgOperatorCartesianSense::mult_M : Unable to sample data: " << 
      cudaGetErrorString(err) << std::endl;
    return -1;
  }

  return 0;
}

int cgOperatorCartesianSense::mult_MH(cuNDArray<float2>* in, cuNDArray<float2>* out, bool accumulate)
{

  if (!(out->dimensions_equal(dimensions_)) ||
      !(in->dimensions_equal(dimensions_out_)) ) {
    std::cerr << "cgOperatorCartesianSense::mult_MH dimensions mismatch" << std::endl;
    return -1;
  }

  if (!accumulate) clear(out);

  return 0;
}

int cgOperatorCartesianSense::mult_MH_M(cuNDArray<float2>* in, cuNDArray<float2>* out, bool accumulate)
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

  if (mult_M(&tmp, out, accumulate) < 0) {
    std::cerr << "cgOperatorCartesianSense::mult_MH_M: Unable to perform mult_M" << std::endl;
    return -2;
  }

  return 0;
}

int cgOperatorCartesianSense::clear(cuNDArray<float2>* in)
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

