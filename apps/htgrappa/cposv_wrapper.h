#ifndef CPOSV_WRAPPER_H
#define CPOSV_WRAPPER_H

#include "cuNDArray.h"

int cposv_wrapper(cuNDArray<float2>* A, cuNDArray<float2>* rhs);


#endif //CPOSV_WRAPPER
