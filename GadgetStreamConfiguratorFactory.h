#ifndef GADGETSTREAMCONFIGURATORFACTOR_H
#define GADGETSTREAMCONFIGURATORFACTOR_H

#include "GadgetStreamConfigurator.h"

class GadgetStreamConfiguratorFactory
{
 public:
  static GadgetStreamConfigurator* 
    CreateConfigurator(GadgetMessageConfigurator c, char* config, GadgetStreamController* controller);
};

#endif //GADGETSTREAMCONFIGURATORFACTOR_H
