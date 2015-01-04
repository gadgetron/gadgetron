#include "ParameterRelayGadget.h"

namespace Gadgetron{
int ParameterRelayGadget
::process(ACE_Message_Block* m)
{
  if (this->next()->putq(m) == -1) {
    m->release();
    GERROR("ParameterRelayGadget::process, passing data on to next gadget");
    return -1;
  }

  return GADGET_OK;
}
GADGET_FACTORY_DECLARE(ParameterRelayGadget)
}


