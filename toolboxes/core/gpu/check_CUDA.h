#pragma once

#include "GadgetronCuException.h"

#include <boost/exception/detail/attribute_noreturn.hpp>

/**
 *  Should never be used in the code, use CHECK_FOR_CUDA_ERROR(); instead
 *  inspired by cutil.h: CUT_CHECK_ERROR
 */
inline void CHECK_FOR_CUDA_ERROR(char const * cur_fun, const char* file, const int line) {
    cudaError_t errorCode = cudaGetLastError();
    if (errorCode != cudaSuccess) {
    	boost::exception_detail::throw_exception_(Gadgetron::cuda_error(errorCode),cur_fun,file,line);
    }
#ifdef DEBUG
    cudaThreadSynchronize();
    errorCode = cudaGetLastError();
    if (errorCode != cudaSuccess) {
    	boost::exception_detail::throw_exception_(Gadgetron::cuda_error(errorCode),cur_fun,file,line);
    }
#endif
}

/**
 *  Checks for CUDA errors and throws an exception if an error was detected.
 */

#define CHECK_FOR_CUDA_ERROR(); CHECK_FOR_CUDA_ERROR(BOOST_CURRENT_FUNCTION,__FILE__,__LINE__);

#define CUDA_CALL(res) {cudaError_t errorCode = res; if (errorCode != cudaSuccess) { boost::exception_detail::throw_exception_(cuda_error(errorCode),BOOST_CURRENT_FUNCTION,__FILE__,__LINE__); } }
