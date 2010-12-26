#include "GadgetStreamConfiguratorFactory.h"

#include "DefaultConfigurator.h"

GadgetStreamConfigurator* 
GadgetStreamConfiguratorFactory
::CreateConfigurator(GadgetMessageConfigurator c, char* config, GadgetStreamController* controller)
{

  ACE_DEBUG( (LM_DEBUG, 
	      ACE_TEXT("ConfiguratorFactory: Creating configurator %s@%s\n"), 
	      c.configurator_name, c.configurator_lib) );

  return new DefaultConfigurator(config, c.configuration_length, controller);
}
