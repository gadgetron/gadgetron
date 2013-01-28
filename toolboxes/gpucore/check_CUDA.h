#pragma once

#include <stdio.h>
#include <boost/exception/detail/attribute_noreturn.hpp>
#include "GadgetronCuException.h"

inline void THROW_ERROR(const char* file, const int line,
                        const char* errorString) {
    const int bLength = 256;
    char buffer[bLength];
    int n = sprintf (buffer,"[file: %s line: %i] CUDA Error: %s\n",
                     file, line,
                     errorString);
    printf( "%s", buffer );
    if (n < 0)
        printf("error when writing error CUDA message\n");
    if (n >= bLength)
        printf("error message buffer was to small\n");
    exit(-1);
}

/**
 *  Should never be used in the code, use CHECK_FOR_CUDA_ERROR(); instead
 *  inspired by cutil.h: CUT_CHECK_ERROR
 */
inline void CHECK_FOR_CUDA_ERROR(char const * cur_fun, const char* file, const int line) {
    cudaError_t errorCode = cudaGetLastError();
    if (errorCode != cudaSuccess) {
    	boost::exception_detail::throw_exception_(cuda_error(errorCode),cur_fun,file,line);
    }
#ifdef DEBUG
    errorCode = cudaThreadSynchronize();
    if (errorCode != cudaSuccess) {
    	boost::exception_detail::throw_exception_(cuda_error(errorCode),cur_fun,file,line);
    }
#endif
}

/**
 *  Checks for CUDA errors and throws an exception if
 *  an error was detected, is only available in debug mode.
 */
//#if DEBUG
#define CHECK_FOR_CUDA_ERROR(); CHECK_FOR_CUDA_ERROR(BOOST_CURRENT_FUNCTION,__FILE__,__LINE__);
//#else
//#define CHECK_FOR_CUDA_ERROR();
//#endif

#define CUDA_CALL(res) {cudaError_t errorCode = res; if (errorCode != cudaSuccess) { boost::exception_detail::throw_exception_(cuda_error(errorCode),BOOST_CURRENT_FUNCTION,__FILE__,__LINE__); } }
