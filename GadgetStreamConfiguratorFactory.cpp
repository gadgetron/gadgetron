#include "ace/OS_NS_string.h"

#include "Gadgetron.h"
#include "GadgetStreamConfiguratorFactory.h"
#include "DefaultConfigurator.h"
#include "GPUCGConfigurator.h"

GadgetStreamConfigurator* 
GadgetStreamConfiguratorFactory
::CreateConfigurator(GadgetMessageConfigurator c, char* config, GadgetStreamController* controller)
{

  ACE_DEBUG( (LM_DEBUG, 
	      ACE_TEXT("ConfiguratorFactory: Creating configurator %s@%s\n"), 
	      c.configurator_name, c.configurator_lib) );

  if (ACE_OS::strcmp(c.configurator_name,"default") == 0 &&
      ACE_OS::strcmp(c.configurator_lib ,"core") == 0) {
    return new DefaultConfigurator(config, c.configuration_length, controller);
  } else if (ACE_OS::strcmp(c.configurator_name,"gpucg") == 0 &&
	     ACE_OS::strcmp(c.configurator_lib ,"gpucg") == 0) {
    return new GPUCGConfigurator(config, c.configuration_length, controller);
  }


  GADGET_DEBUG2("Unrecognized configurator %s@%s, returning deault configurator\n",
		c.configurator_name, c.configurator_lib);

  return new DefaultConfigurator(config, c.configuration_length, controller);
}
