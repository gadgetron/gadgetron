#pragma once

#include <stdio.h>

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
inline void CHECK_FOR_CUDA_ERROR(const char* file, const int line) {
    cudaError_t errorCode = cudaGetLastError();
    if (errorCode != cudaSuccess) {
        const char* errorString = cudaGetErrorString(errorCode);
        THROW_ERROR(file, line, errorString);
    }
#ifdef DEBUG
    errorCode = cudaThreadSynchronize();
    if (errorCode != cudaSuccess) {
        const char* errorString = cudaGetErrorString(errorCode);
        THROW_ERROR(file, line, errorString);
    }
#endif
}

/**
 *  Checks for CUDA errors and throws an exception if
 *  an error was detected, is only available in debug mode.
 */
//#if DEBUG
#define CHECK_FOR_CUDA_ERROR(); CHECK_FOR_CUDA_ERROR(__FILE__,__LINE__);
//#else
//#define CHECK_FOR_CUDA_ERROR();
//#endif
