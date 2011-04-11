#include "GPUTimer.h"

#include <iostream>

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
  double time = cutGetTimerValue( timer_ ); 
  std::cout << name_ << ": " << time << " ms" << std::endl;
}
