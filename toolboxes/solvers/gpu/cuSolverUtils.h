#pragma once
#include "complext.h"
#include "hoCuNDArray.h"
#include "cuNDArray.h"
#include "gpusolvers_export.h"

namespace Gadgetron{

template<class T> void EXPORTGPUSOLVERS solver_non_negativity_filter(cuNDArray<T>* x , cuNDArray<T>* g);

}
