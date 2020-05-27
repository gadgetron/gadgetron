/**
 * \file cuSDC_kernel.h
 * \brief Kernel for sampling density compensation (CUDA specialization).
 */

#pragma once

#include "SDC_kernel.h"

#include "cuNDArray.h"

namespace Gadgetron
{
    template<class REAL, unsigned int D>
    using cuSDC_kernel = SDC_kernel<cuNDArray, REAL, D>;
}