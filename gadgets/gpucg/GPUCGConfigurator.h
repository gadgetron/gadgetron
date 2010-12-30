#ifndef GPUCGCONFIGURATOR_H
#define GPUCGCONFIGURATOR_H

#include "GadgetStreamConfigurator.h"

class GPUCGConfigurator : public GadgetStreamConfigurator
{
 public:
  GPUCGConfigurator(char* config, ACE_UINT16 config_len, GadgetStreamController* controller);
  virtual int ConfigureStream(ACE_Stream<ACE_MT_SYNCH>* stream);
};

#endif //GADGETDEFAULTCONFIGURATOR_H
