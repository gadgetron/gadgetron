#include "GadgetronInfoGadget.h"
#include <string>

namespace Gadgetron{
  
  int GadgetronInfoGadget::process(ACE_Message_Block* m)
  {
    if (this->next()->putq(m) == -1) {
      m->release();
      GERROR("GadgetroninfoGadget::process, passing data on to next gadget");
      return -1;
    }
    
    return GADGET_OK;
  }


  int GadgetronInfoGadget::close(unsigned long flags)
  {
    if ( Gadget::close(flags) != GADGET_OK ) return GADGET_FAIL;

    GadgetContainerMessage<std::string>* m1 = new GadgetContainerMessage<std::string>();

    GadgetContainerMessage<GadgetMessageIdentifier>* mb = new GadgetContainerMessage<GadgetMessageIdentifier>();
    mb->getObjectPtr()->id = GADGET_MESSAGE_TEXT;
    mb->cont(m1);
    
    int ret =  this->controller_->output_ready(mb);
    if ( (ret < 0) ) {
      mb->release();
      GDEBUG("Failed to return massage to controller\n");
      return GADGET_FAIL;
    }

    return GADGET_OK;
  }
  
  GADGET_FACTORY_DECLARE(GadgetronInfoGadget)
}
