/**
 * \file hoSDC_kernel.h
 * \brief Kernel for sampling density compensation (host specialization).
 */

#pragma once

#include "SDC_kernel.h"

#include "hoNDArray.h"

namespace Gadgetron
{   
    template<class REAL, unsigned int D>
    using hoSDC_kernel = SDC_kernel<hoNDArray, REAL, D>;
}