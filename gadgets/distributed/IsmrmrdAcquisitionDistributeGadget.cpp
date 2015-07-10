#include "IsmrmrdAcquisitionDistributeGadget.h"
#include "GadgetMRIHeaders.h"

namespace Gadgetron{

  int IsmrmrdAcquisitionDistributeGadget::node_index(ACE_Message_Block* m)
  {
    auto h = AsContainerMessage<ISMRMRD::AcquisitionHeader>(m);

    if (!h) return GADGET_FAIL;
    
    return h->getObjectPtr()->idx.slice;
  }

  int IsmrmrdAcquisitionDistributeGadget::message_id(ACE_Message_Block* m)
  {
    auto h = AsContainerMessage<ISMRMRD::AcquisitionHeader>(m);
    if (!h) return 0;

    return GADGET_MESSAGE_ISMRMRD_ACQUISITION;
  }

  GADGET_FACTORY_DECLARE(IsmrmrdAcquisitionDistributeGadget)

}


