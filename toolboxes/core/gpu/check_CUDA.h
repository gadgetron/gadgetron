/** \file check_CUDA.h
    \brief Macroes to check whether GPU-based code has caused any errors, and if so, throw a runtime exception accordingly.
*/

#pragma once

#include "GadgetronCuException.h"

namespace Gadgetron {

  /**
   *  Should never be used in the code, use CHECK_FOR_CUDA_ERROR(); instead
   *  inspired by cutil.h: CUT_CHECK_ERROR
   */
  inline void CHECK_FOR_CUDA_ERROR(char const * cur_fun, const char* file, const int line) {
    cudaError_t errorCode = cudaGetLastError();
    if (errorCode != cudaSuccess) {
      throw cuda_error(errorCode);
    }
#ifdef DEBUG
    cudaThreadSynchronize();
    errorCode = cudaGetLastError();
    if (errorCode != cudaSuccess) {
      throw cuda_error(errorCode);
    }
#endif
  }
}

/**
 *  Checks for CUDA errors and throws an exception if an error was detected.
 */
#define CHECK_FOR_CUDA_ERROR(); CHECK_FOR_CUDA_ERROR(BOOST_CURRENT_FUNCTION,__FILE__,__LINE__);

/**
 *  Call "res", checks for CUDA errors and throws an exception if an error was detected.
 */
#define CUDA_CALL(res) {cudaError_t errorCode = res; if (errorCode != cudaSuccess) { throw cuda_error(errorCode); }}
