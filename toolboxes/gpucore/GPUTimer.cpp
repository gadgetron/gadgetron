#include "GPUTimer.h"

#include <iostream>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cutil.h>

GPUTimer::GPUTimer(const char* name)
  : name_(name)
{
  cutCreateTimer(&timer_);
  cutResetTimer( timer_ ); cutStartTimer( timer_ );
}

GPUTimer::GPUTimer()
{
  GPUTimer("GPUTimer");
}

GPUTimer::~GPUTimer()
{
  cudaThreadSynchronize();
  double time = cutGetTimerValue( timer_ ); 
  std::cout << name_ << ": " << time << " ms" << std::endl;
}
