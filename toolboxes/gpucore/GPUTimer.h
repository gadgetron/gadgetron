#ifndef __GPUTIMER_H
#define __GPUTIMER_H

#pragma once
#include "gadgetron_export.h"

#include <string>
#include <cutil.h>

class GPUTimer
{
 public:

  GPUTimer() { GPUTimer("GPUTimer"); }

  GPUTimer(const char* name) : name_(name) {
  cutCreateTimer(&timer_); cutResetTimer( timer_ ); cutStartTimer( timer_ );
}


  virtual ~GPUTimer() {
    cudaThreadSynchronize();
    double time = cutGetTimerValue( timer_ ); 
    std::cout << name_ << ": " << time << " ms" << std::endl;
  }

 protected:
  unsigned int timer_;
  std::string name_;
};

#endif //__GPUTIMER_H
