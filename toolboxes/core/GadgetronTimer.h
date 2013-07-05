/** \file GadgetronTimer.h
    \brief Generic timer class to measure runtime performance.
*/

#ifndef __GADGETRONTIMER_H
#define __GADGETRONTIMER_H

#pragma once

#ifdef WIN32 
#include <windows.h>
#else 
#include <sys/time.h>
#endif

#include <iostream>
#include <string>

namespace Gadgetron{

  class GadgetronTimer
  {
  public:

    GadgetronTimer() { GadgetronTimer("GPUTimer"); }

    GadgetronTimer(const char* name) : name_(name) 
    {
      pre();
#ifdef WIN32
      QueryPerformanceFrequency(&frequency_);
      QueryPerformanceCounter(&start_);
#else
      gettimeofday(&start_, NULL);
#endif
    }
    
    virtual ~GadgetronTimer() 
    {
      double time_in_us = 0.0;
      post();
#ifdef WIN32
      QueryPerformanceCounter(&end_);
      time_in_us = (end_.QuadPart * (1.0e6/ frequency_.QuadPart)) - start_.QuadPart * (1.0e6 / frequency_.QuadPart);
#else
      gettimeofday(&end_, NULL);
      time_in_us = ((end_.tv_sec * 1e6) + end_.tv_usec) - ((start_.tv_sec * 1e6) + start_.tv_usec);
#endif
      std::cout << name_ << ": " << time_in_us/1000.0 << " ms" << std::endl; std::cout.flush();
    }
    
    virtual void pre() {}
    virtual void post() {}
    
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
}

#endif //__GADGETRONTIMER_H
