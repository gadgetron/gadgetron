#ifndef GADGETDEFAULTCONFIGURATOR_H
#define GADGETDEFAULTCONFIGURATOR_H

#include "GadgetStreamConfigurator.h"

class DefaultConfigurator : public GadgetStreamConfigurator
{
 public:
  DefaultConfigurator(char* config, ACE_UINT16 config_len);
  virtual int ConfigureStream(ACE_Stream<ACE_MT_SYNCH>* stream);
};

#endif //GADGETDEFAULTCONFIGURATOR_H
