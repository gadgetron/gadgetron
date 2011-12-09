#ifndef __GPUTIMER_H
#define __GPUTIMER_H

#pragma once

#include <cuda_runtime_api.h>
#include "../gadgettools/GadgetronTimer.h"

class GPUTimer : public GadgetronTimer
{
public:
	GPUTimer()
		: GadgetronTimer()
	{

	}

	GPUTimer(const char* name)
		: GadgetronTimer(name)
	{

	}

	virtual void pre() {
	    cudaThreadSynchronize();
	}

	virtual void post() {
		cudaThreadSynchronize();
	}
};

#endif //__GPUTIMER_H
