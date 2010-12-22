#ifndef GADGETSTREAMCONFIGURATOR_H
#define GADGETSTREAMCONFIGURATOR_H

#include "ace/Stream.h"

#include "gadgetheaders.h"

class GadgetStreamConfigurator
{
 public:
  GadgetStreamConfigurator(char* config, ACE_UINT16 config_len);

  virtual ~GadgetStreamConfigurator();
  virtual int ConfigureStream(ACE_Stream<ACE_MT_SYNCH>* stream) = 0;
 protected:
  ACE_TCHAR* config_;
  ACE_UINT16 config_length_;

};

GadgetStreamConfigurator* 
CreateGadgetStreamConfigurator(GadgetMessageConfigurator, char* config);

#endif //GADGETSTREAMCONFIGURATOR_H
