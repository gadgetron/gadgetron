#ifndef GADGETHEADERS_H
#define GADGETHEADERS_H

#include <ace/Basic_Types.h>

enum GadgetMessageID {
  GADGET_MESSAGE_ID_MIN = 0,
  GADGET_MESSAGE_ACQUISITION,
  GADGET_MESSAGE_CONFIGURATION,
  GADGET_MESSAGE_NEW_MEASUREMENT,
  GADGET_MESSAGE_END_OF_SCAN,
  GADGET_MESSAGE_IMAGE,
  GADGET_MESSAGE_EMPTY,
  GADGET_MESSAGE_ID_MAX
};
  
struct GadgetMessageIdentifier
{
  ACE_UINT16 id;
};

struct GadgetMessageConfigurator
{
  ACE_UINT16 config_id;
  char configuration_info[1024];
};


#endif  //GADGETHEADERS_H
