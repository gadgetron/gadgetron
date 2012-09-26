#include "GadgetIsmrmrdReadWrite.h"
#include "GadgetMessageInterface.h"
#include "AcquisitionFinishGadget.h"
#include "GadgetStreamController.h"

int AcquisitionFinishGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
				 GadgetContainerMessage< NDArray< std::complex<float> > >* m2)
{
  if (!controller_) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Cannot return result to controller, no controller set")) );
    return -1;
  }

  GadgetContainerMessage<GadgetMessageIdentifier>* mb =
    new GadgetContainerMessage<GadgetMessageIdentifier>();

  mb->getObjectPtr()->id = GADGET_MESSAGE_ACQUISITION;

  mb->cont(m1);

  return controller_->output_ready(mb);

}

GADGET_FACTORY_DECLARE(AcquisitionFinishGadget)
