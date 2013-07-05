#pragma once

#include <cuda_runtime_api.h>
#include <stdexcept>

namespace Gadgetron{
  
  class cuda_error : public std::runtime_error
  {
  public:
    cuda_error(std::string msg) : std::runtime_error(msg) {}
    cuda_error(cudaError_t errN) : std::runtime_error(cudaGetErrorString(errN)) {
    }
  };
}
