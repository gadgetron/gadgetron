#include "GadgetronInfoGadget.h"
#include "gadgetron_system_info.h"

#include <string>
#include <sstream>

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
    if ( BasicPropertyGadget::close(flags) != GADGET_OK ) return GADGET_FAIL;

    if (!info_sent_) {
      info_sent_ = true;
      
      std::stringstream ss;
      print_system_information(ss);
      ACE_Message_Block* m1 = new ACE_Message_Block(ss.str().size());
      memcpy(m1->wr_ptr(),ss.str().c_str(),ss.str().size());
      
      GadgetContainerMessage<GadgetMessageIdentifier>* mb = new GadgetContainerMessage<GadgetMessageIdentifier>();
      mb->getObjectPtr()->id = GADGET_MESSAGE_TEXT;
      mb->cont(m1);
      
      int ret =  this->controller_->output_ready(mb);
      if ( (ret < 0) ) {
	mb->release();
	GDEBUG("Failed to return massage to controller\n");
	return GADGET_FAIL;
      }
    }
    return GADGET_OK;
  }
  
  GADGET_FACTORY_DECLARE(GadgetronInfoGadget)
}
