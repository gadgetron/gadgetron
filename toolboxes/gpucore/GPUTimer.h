#ifndef __GPUTIMER_H
#define __GPUTIMER_H

#pragma once
#include "gadgetron_export.h"

#include <string>

class EXPORTGPUCORE GPUTimer
{
 public:
  GPUTimer(const char* name);
  GPUTimer();
  virtual ~GPUTimer();

 protected:
  unsigned int timer_;
  std::string name_;
};

#endif //__GPUTIMER_H
