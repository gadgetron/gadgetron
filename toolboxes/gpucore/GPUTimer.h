#ifndef __GPUTIMER_H
#define __GPUTIMER_H

#pragma once

#ifdef WIN32 
#include <windows.h>
#else 
#include <sys/time.h>
#endif

#include <iostream>
#include <cuda_runtime_api.h>

class GPUTimer
{
 public:

  GPUTimer() { GPUTimer("GPUTimer"); }

  GPUTimer(const char* name) : name_(name) {
    cudaThreadSynchronize(); //Should we keep this here?
#ifdef WIN32
    QueryPerformanceFrequency(&frequency_);
    QueryPerformanceCounter(&start_);
#else
    gettimeofday(&start_, NULL);
#endif
  }


  virtual ~GPUTimer() {
    cudaThreadSynchronize(); //Should we keep this here, if it is here, it should also be in the constructor
    double time_in_us = 0.0;
#ifdef WIN32
    QueryPerformanceCounter(&end_);
    endTimeInMicroSec = (end_.QuadPart * (1.0e6/ frequency_.QuadPart)) - 
      start_.QuadPart * (1.0e6 / frequency_.QuadPart);
#else
    gettimeofday(&end_, NULL);
    time_in_us = ((end_.tv_sec * 1e6) + end_.tv_usec) - ((start_.tv_sec * 1e6) + start_.tv_usec);
#endif
    std::cout << name_ << ": " << time_in_us/1000.0 << " ms" << std::endl;
  }

 protected:
#ifdef WIN32
    LARGE_INTEGER frequency_;
    LARGE_INTEGER start_;
    LARGE_INTEGER end_;
#else
    timeval start_;
    timeval end_;
#endif

  std::string name_;
};

#endif //__GPUTIMER_H
