#include "IsmrmrdAcquisitionDistributeGadget.h"

namespace Gadgetron{

  int IsmrmrdAcquisitionDistributeGadget::node_index(ACE_Message_Block* m)
  {
    return 0;
  }
  
  GADGET_FACTORY_DECLARE(IsmrmrdAcquisitionDistributeGadget)

}


