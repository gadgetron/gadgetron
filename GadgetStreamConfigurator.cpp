#include "GadgetStreamConfigurator.h"

GadgetStreamConfigurator
::GadgetStreamConfigurator(char* config, ACE_UINT16 config_len)
  : config_(0)
  , config_length_(0)
{
  config_length_ = config_len;
  config_ = new ACE_TCHAR[config_length_];
  if (!config_) {
    ACE_DEBUG( (LM_ERROR, 
		ACE_TEXT("Unable to allocate memory for config info\n")) );
    config_length_ = 0;
  }
  
  ACE_OS::strncpy(config_,config,config_length_);

}

GadgetStreamConfigurator::~GadgetStreamConfigurator()
{
  if (config_) delete [] config_;
}
