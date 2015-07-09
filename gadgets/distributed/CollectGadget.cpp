#include "CollectGadget.h"

namespace Gadgetron{

  int CollectGadget::process(ACE_Message_Block* m)
  {
    //It is enough to put the first one, since they are linked
    if (this->next()->putq(m) == -1) {
      m->release();
      GERROR("CollectGadget::process, passing data on to next gadget");
      return -1;
    }
    
    return 0;
  }

  GADGET_FACTORY_DECLARE(CollectGadget)
}


