#include "cuNDFFT.h"

#include <cufft.h>
#include <cublas.h>

int cuNDFFT::fft_int(cuNDArray< cuFloatComplex >* input, std::vector<unsigned int> dims_to_transform, 
		     int direction, bool do_scale)
{
  std::vector<unsigned int> new_dim_order;
  std::vector<unsigned int> reverse_dim_order;
  std::vector<int> dims;
  std::vector<unsigned int> dim_count(input->get_number_of_dimensions(),0);
  
  unsigned int array_ndim = input->get_number_of_dimensions();
  std::vector<unsigned int> array_dims = input->get_dimensions();
  
  dims = std::vector<int>(dims_to_transform.size(),0);
  for (unsigned int i = 0; i < dims_to_transform.size(); i++) {
    if (dims_to_transform[i] >= array_ndim) {
      std::cerr << "cuNDFFT::fft Invalid dimensions specified for transform " << dims_to_transform[i] << "max " << array_ndim << std::endl;
      return -1;
    }
    if (dim_count[dims_to_transform[i]] > 0) {
      std::cerr << "cuNDFFT::fft Invalid dimensions (duplicates) specified for transform" << std::endl;
      return -1;
    }
    dim_count[dims_to_transform[i]]++;
    dims[dims_to_transform.size()-1-i] = array_dims[dims_to_transform[i]];
  }
  
  new_dim_order = dims_to_transform;
  for (unsigned int i = 0; i < array_ndim; i++) {
    if (!dim_count[i]) new_dim_order.push_back(i);
  }

  reverse_dim_order = std::vector<unsigned int>(array_ndim,0);
  for (unsigned int i = 0; i < array_ndim; i++) {
    reverse_dim_order[new_dim_order[i]] = i;
  }

  int ndim = dims.size();
  int batches = 0;
  int elements_in_ft = 1;
  for (unsigned int i = 0; i < dims.size(); i++) elements_in_ft *= dims[i];
  batches = input->get_number_of_elements() / elements_in_ft;


  cufftHandle plan;
  cufftResult ftres;

  ftres = cufftPlanMany(&plan,ndim,&dims[0],&dims[0],1,elements_in_ft,&dims[0],1,elements_in_ft,CUFFT_C2C,batches);
  if (ftres != CUFFT_SUCCESS) {
    std::cerr << "cuNDFFT FFT plan failed: " << ftres << std::endl;
    return -1;
  }

  //IFFTSHIFT
  if (input->permute(new_dim_order,0,-1) < 0) {
    std::cerr << "cuNDFFT error permuting before FFT" << std::endl;
    return -1;
  }

  if (cufftExecC2C(plan, input->get_data_ptr(), input->get_data_ptr(), direction) != CUFFT_SUCCESS) {
    std::cerr << "cuNDFFT FFT execute failed" << std::endl;
    return -1;

  }

  if (do_scale) {
    cuFloatComplex scale;
    scale.x = 1.0f/elements_in_ft;
    scale.y = 0.0f;
    cublasCscal(input->get_number_of_elements(), scale, input->get_data_ptr(), 1);
  }

  //FFTSHIFT 
  if (input->permute(reverse_dim_order,0,1) < 0) {
    std::cerr << "cuNDFFT error permuting after FFT" << std::endl;
    return -1;
  }
  
  return 0;
}

int cuNDFFT::fft(cuNDArray< cuFloatComplex >* input, std::vector<unsigned int> dims_to_transform)
{
  return fft_int(input,dims_to_transform, CUFFT_FORWARD, false);
}

int cuNDFFT::ifft(cuNDArray< cuFloatComplex >* input, std::vector<unsigned int> dims_to_transform, bool do_scale)
{
  return fft_int(input,dims_to_transform, CUFFT_INVERSE, do_scale);
}


int cuNDFFT::fft(cuNDArray< cuFloatComplex >* input, unsigned int dim_to_transform)
{
  std::vector<unsigned int> dims(1,dim_to_transform);
  return fft_int(input,dims, CUFFT_FORWARD, false);
}
  
int cuNDFFT::ifft(cuNDArray< cuFloatComplex >* input, unsigned int dim_to_transform, bool do_scale)
{
  std::vector<unsigned int> dims(1,dim_to_transform);
  return fft_int(input,dims, CUFFT_INVERSE, do_scale);
}

int cuNDFFT::fft(cuNDArray< cuFloatComplex >* input)
{
  std::vector<unsigned int> dims(input->get_number_of_dimensions(),0);
  for (unsigned int i = 0; i < dims.size(); i++) dims[i] = i;

  return fft_int(input,dims, CUFFT_FORWARD, false);
}

int cuNDFFT::ifft(cuNDArray< cuFloatComplex >* input, bool do_scale)
{
  std::vector<unsigned int> dims(input->get_number_of_dimensions(),0);
  for (unsigned int i = 0; i < dims.size(); i++) dims[i] = i;
  return fft_int(input,dims, CUFFT_INVERSE, do_scale);
}

