#pragma once

#include "GadgetronException.h"
#include <cuda_runtime_api.h>

namespace Gadgetron{
  
  class cuda_error : public runtime_error
  {
  public:
    cuda_error(std::string msg) : runtime_error(msg) {}
    cuda_error(cudaError_t errN) : runtime_error(cudaGetErrorString(errN)) {
    }
  };
}
