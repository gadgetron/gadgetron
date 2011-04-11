#ifndef __GPUTIMER_H
#define __GPUTIMER_H

#include <string>

class GPUTimer
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
