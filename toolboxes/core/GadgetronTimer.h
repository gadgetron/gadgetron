/** \file GadgetronTimer.h
    \brief Generic timer class to measure runtime performance.
*/

#ifndef __GADGETRONTIMER_H
#define __GADGETRONTIMER_H

#pragma once

#include <sys/time.h>

#include <string>
#include "log.h"

namespace Gadgetron{

  class GadgetronTimer
  {
  public:

    GadgetronTimer() : name_("GPUTimer"), timing_in_destruction_(true)
    {
        pre();
        start();
    }

    GadgetronTimer(bool timing) : name_("GPUTimer"), timing_in_destruction_(timing)
    {
        if ( timing_in_destruction_ )
        {
            pre();
            start();
        }
    }

    GadgetronTimer(const char* name, bool timing=true) : name_(name), timing_in_destruction_(timing) 
    {
        if ( timing_in_destruction_ )
        {
            pre();
            start();
        }
    }

    virtual ~GadgetronTimer() 
    {
        if ( timing_in_destruction_ )
        {
            post();
            stop();
        }
    }

    virtual void pre() {}
    virtual void post() {}

    virtual void start()
    {
        gettimeofday(&start_, NULL);
    }

    void start(const char* name)
    {
        name_ = name;
        start();
    }

    virtual double stop()
    {
        double time_in_us = 0.0;
        gettimeofday(&end_, NULL);
        time_in_us = ((end_.tv_sec * 1e6) + end_.tv_usec) - ((start_.tv_sec * 1e6) + start_.tv_usec);
        GDEBUG("%s:%f ms\n", name_.c_str(), time_in_us/1000.0);
        return time_in_us;
    }

    void set_timing_in_destruction(bool timing) { timing_in_destruction_ = timing; }

  protected:

    timeval start_;
    timeval end_;

    std::string name_;

    bool timing_in_destruction_;
  };
}

#endif //__GADGETRONTIMER_H
