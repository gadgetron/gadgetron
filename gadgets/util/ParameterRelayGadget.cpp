#include "ParameterRelayGadget.h"
#include "Gadgetron.h"
namespace Gadgetron{
int ParameterRelayGadget
::process(ACE_Message_Block* m)
{
  if (this->next()->putq(m) == -1) {
    m->release();
    ACE_ERROR_RETURN( (LM_ERROR,
		       ACE_TEXT("%p\n"),
		       ACE_TEXT("ParameterRelayGadget::process, passing data on to next gadget")),
		      -1);
  }

  return GADGET_OK;
}
GADGET_FACTORY_DECLARE(ParameterRelayGadget)
}


