#include "CollectGadget.h"
#include <ismrmrd/ismrmrd.h>
#include "GadgetMRIHeaders.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "GadgetStreamInterface.h"

namespace Gadgetron{

  int CollectGadget::message_id(ACE_Message_Block* m)
  {
    if (AsContainerMessage<ISMRMRD::AcquisitionHeader>(m)) {
      return GADGET_MESSAGE_ISMRMRD_ACQUISITION;
    } else if (AsContainerMessage<ISMRMRD::ImageHeader>(m)) {
      return GADGET_MESSAGE_ISMRMRD_IMAGE;
    } else {
      return 0;
    }
  }

  CollectGadget::CollectGadget()
  {
  }

  CollectGadget::~CollectGadget()
  {
  }

  int CollectGadget::process(ACE_Message_Block* m)
  {
    if (pass_through_mode.value()) {
      //It is enough to put the first one, since they are linked
      if (this->next()->putq(m) == -1) {
        m->release();
        GERROR("CollectGadget::process, passing data on to next gadget");
        return -1;
      }
    } else {
      if (!this->controller_)
      {
        GERROR("Cannot return result to controller, no controller set");
        return -1;
      }

      GadgetContainerMessage<GadgetMessageIdentifier>* mb = new GadgetContainerMessage<GadgetMessageIdentifier>();

      mb->getObjectPtr()->id = message_id(m);
      mb->cont(m);

      int ret = this->controller_->output_ready(mb);

      if ((ret < 0))
      {
        GERROR("Failed to return massage to controller\n");
        return GADGET_FAIL;
      }

      return GADGET_OK;
    }

    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(CollectGadget)
}
