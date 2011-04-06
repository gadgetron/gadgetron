#include "cuFFT.h"

int cuFFT::fft(cuNDArray< cuFloatComplex >* input, unsigned int dim_to_transform)
{

  /*
  if (!is_array_valid(input)) {
    return -1;
  }
  */

  /*
  K2I( input->get_data_ptr(), 
       uintvec_to_uint4(input->get_dimensions()), 
       dim_to_transform);
  */

  return 0;
}
  
int cuFFT::ifft(cuNDArray< cuFloatComplex >* input, unsigned int dim_to_transform)
{
  /*
  if (!is_array_valid(input)) {
    return -1;
  }
  */

  /*

  K2I( input->get_data_ptr(), 
       uintvec_to_uint4(input->get_dimensions()), 
       dim_to_transform);
  */

  return 0;
}

int cuFFT::fft(cuNDArray< cuFloatComplex >* input)
{
  

  return 0;
}

int cuFFT::ifft(cuNDArray< cuFloatComplex >* input)
{
  return 0;
}

bool is_array_valid(cuNDArray< cuFloatComplex >* in)
{
  if (in->get_number_of_dimensions() > 4) {
    std::cerr << "cuFFT: arrays with more dimensions than 4 are not supported at the moment" 
	      << std::endl;

    return false;
  }

  return true;
}

uint2 cuFFT::uintvec_to_uint2(std::vector<unsigned int>& vec)
{
  uint2 ret;

  if (vec.size() < 1) 
    ret.x = 1;
  else 
    ret.x = vec[0];

  if (vec.size() < 2) 
    ret.y = 1;
  else 
    ret.y = vec[1];

  return ret;
}

uint3 cuFFT::uintvec_to_uint3(std::vector<unsigned int>& vec)
{
  uint3 ret;

  if (vec.size() < 1) 
    ret.x = 1;
  else 
    ret.x = vec[0];

  if (vec.size() < 2) 
    ret.y = 1;
  else 
    ret.y = vec[1];

  if (vec.size() < 3) 
    ret.z = 1;
  else 
    ret.z = vec[2];

  return ret;
}

uint4 cuFFT::uintvec_to_uint4(std::vector<unsigned int>& vec)
{
  uint4 ret;

  if (vec.size() < 1) 
    ret.x = 1;
  else 
    ret.x = vec[0];

  if (vec.size() < 2) 
    ret.y = 1;
  else 
    ret.y = vec[1];

  if (vec.size() < 3) 
    ret.z = 1;
  else 
    ret.z = vec[2];

  if (vec.size() < 4) 
    ret.w = 1;
  else 
    ret.w = vec[3];

  return ret;
}
