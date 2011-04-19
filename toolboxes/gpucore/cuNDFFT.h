#ifndef CUNDFFT_H
#define CUNDFFT_H

#include "cuda.h"
#include "cuComplex.h"
#include "cuNDArray.h"

class cuNDFFT
{

 public:
  cuNDFFT() {}

  int fft(cuNDArray< cuFloatComplex >* input, std::vector<unsigned int> dims_to_transform);
  int ifft(cuNDArray< cuFloatComplex >* input, std::vector<unsigned int> dims_to_transform, bool do_scale = true);

  int fft(cuNDArray< cuFloatComplex >* input, unsigned int dim_to_transform);
  int ifft(cuNDArray< cuFloatComplex >* input, unsigned int dim_to_transform, bool do_scale = true);

  int fft(cuNDArray< cuFloatComplex >* input);
  int ifft(cuNDArray< cuFloatComplex >* input, bool do_scale = true);

 protected:
  int fft_int(cuNDArray< cuFloatComplex >* input, std::vector<unsigned int> dims_to_transform, int direction, bool do_scale = true);
};

#endif
