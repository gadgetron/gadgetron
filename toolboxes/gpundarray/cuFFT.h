#ifndef CUFFT_H
#define CUFFT_H

#include "cuda.h"
#include "cuComplex.h"

#include "cuNDArray.h"
#include "FFT.hcu"

class cuFFT
{

 public:
  cuFFT() {}

  int fft(cuNDArray< cuFloatComplex >* input, unsigned int dim_to_transform);
  int ifft(cuNDArray< cuFloatComplex >* input, unsigned int dim_to_transform);
  int fft(cuNDArray< cuFloatComplex >* input);
  int ifft(cuNDArray< cuFloatComplex >* input);

 protected:
  bool is_array_valid(cuNDArray< cuFloatComplex >* in);
  uint2 uintvec_to_uint2(std::vector<unsigned int>& vec);
  uint3 uintvec_to_uint3(std::vector<unsigned int>& vec);
  uint4 uintvec_to_uint4(std::vector<unsigned int>& vec);

};

#endif
