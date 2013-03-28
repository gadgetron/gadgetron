#pragma once

#include "cuNDArray.h"

namespace Gadgetron{


template<class T> void impulse_kernel(cuNDArray<T>* x, cuNDArray<T>* histogram, T&, T& );


}
