#ifndef GADGETRONINFOGADGET_H
#define GADGETRONINFOGADGET_H

#include "Gadget.h"
#include "gadgetron_util_gadgets_export.h"
#include "GadgetMessageInterface.h"
#include "GadgetStreamController.h"

namespace Gadgetron{

  class EXPORTUTILGADGETS GadgetronInfoGadget : public BasicPropertyGadget
    {
    public:
      GadgetronInfoGadget()
	: BasicPropertyGadget()
	{
	  info_sent_ = false;
	}
      GADGET_DECLARE(GadgetronInfoGadget);
      
    protected:
      virtual int process(ACE_Message_Block* m);
      virtual int close(unsigned long flags);
      bool info_sent_;
    };
}
#endif //GADGETRONINFOGADGET_H
