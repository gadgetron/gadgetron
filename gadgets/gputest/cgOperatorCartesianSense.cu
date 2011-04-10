#include "cgOperatorCartesianSense.h"

#include "cuNDFFT.h"

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

__global__ void insert_samples_kernel(float2* in, float2* out, 
				      unsigned int* idx, 
				      unsigned long image_elements,
				      unsigned long int samples,
				      unsigned int coils)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < samples) {
    for (unsigned int i = 0; i < coils; i++) {
      out[idx[idx_in] + i*image_elements].x += in[idx_in + i*samples].x;
      out[idx[idx_in] + i*image_elements].y += in[idx_in + i*samples].y;
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

    return -1;
  }

  cuNDArray<float2> tmp;
  std::vector<unsigned int> full_dimensions = dimensions_;
  full_dimensions.push_back(coils_);

  if (!tmp.create(full_dimensions)) {
    std::cerr << "cgOperatorCartesianSense::mult_M unable to allocate temp array" << std::endl;
    return -1;    
  }

  if (mult_csm(in,&tmp) < 0) {
    std::cerr << "cgOperatorCartesianSense::mult_M : Unable to multiply with coil sensitivities" << std::endl;
    return -1;
  }

  cuNDFFT ft;
  std::vector<unsigned int> ft_dims;
  for (unsigned int i = 0; i < dimensions_.size(); i++) {
    ft_dims.push_back(i);
  }

  ft.fft(&tmp, ft_dims);

  if (!accumulate) clear(out);

  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)idx_->get_number_of_elements()/blockDim.x), 1, 1 );
  sample_array_kernel<<< gridDim, blockDim >>>( tmp.get_data_ptr(), out->get_data_ptr(), idx_->get_data_ptr(),
						in->get_number_of_elements(), idx_->get_number_of_elements(), coils_);
  cudaError_t err = cudaGetLastError();
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

  std::vector<unsigned int> tmp_dimensions = dimensions_;
  tmp_dimensions.push_back(coils_);

  cuNDArray<float2> tmp;
  if (!tmp.create(tmp_dimensions)) {
    std::cerr << "cgOperatorCartesianSense::mult_MH: Unable to create temp storage" << std::endl;
    return -1;
  }

  clear(&tmp);

  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)idx_->get_number_of_elements()/blockDim.x), 1, 1 );
  insert_samples_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), tmp.get_data_ptr(),
						  idx_->get_data_ptr(),out->get_number_of_elements(),
						  idx_->get_number_of_elements(), coils_);
  
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "cgOperatorCartesianSense::mult_EM : Unable to insert samples into array: " << 
      cudaGetErrorString(err) << std::endl;
    return -1;
  }

  cuNDFFT ft;
  std::vector<unsigned int> ft_dims;
  for (unsigned int i = 0; i < dimensions_.size(); i++) {
    ft_dims.push_back(i);
  }

  ft.ifft(&tmp, ft_dims);

  if (!accumulate) clear(out);
  
  if (mult_csm_conj_sum(&tmp,out) < 0) {
    std::cerr << "cgOperatorCartesianSense::mult_MH: Unable to multiply with conjugate of sensitivity maps and sum" << std::endl;
    return -1;
 
  }

  return 0;
}





