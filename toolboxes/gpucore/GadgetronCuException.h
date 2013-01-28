#pragma once

#include "GadgetronException.h"
class cuda_error : public gt_runtime_error
{
public:
	cuda_error(std::string msg) : gt_runtime_error(msg) {}
	cuda_error(cudaError_t errN) : gt_runtime_error(cudaGetErrorString(errN)) {

	}
};
