#include "DefaultConfigurator.h"

DefaultConfigurator::DefaultConfigurator(char* config, ACE_UINT16 config_len)
  : GadgetStreamConfigurator(config,config_len)
{


}

int DefaultConfigurator::ConfigureStream(ACE_Stream<ACE_MT_SYNCH>* stream)
{
  return 0;
}
