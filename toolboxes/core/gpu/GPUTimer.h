/** file GPUTimer.h
    Utility to measure Cuda performance. 
*/

#ifndef __GPUTIMER_H
#define __GPUTIMER_H

#pragma once

#include "GadgetronTimer.h"
#include <cuda_runtime_api.h>

namespace Gadgetron{

    class GPUTimer : public GadgetronTimer
    {
    public:
        GPUTimer() : GadgetronTimer() {}
        GPUTimer(const char* name) : GadgetronTimer(name) {}

        virtual void pre() {
            cudaDeviceSynchronize();
        }

        virtual void post() {
        		cudaDeviceSynchronize();
        }
    };
}
#endif //__GPUTIMER_H
