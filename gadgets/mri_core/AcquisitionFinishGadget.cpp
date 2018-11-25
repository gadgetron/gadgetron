#include "GadgetIsmrmrdReadWrite.h"
#include "GadgetMessageInterface.h"
#include "AcquisitionFinishGadget.h"

using namespace Gadgetron;

int AcquisitionFinishGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
				 GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
	/*
  if (!controller_) {
    GERROR("Cannot return result to controller, no controller set\n");
    return -1;
  }

  GadgetContainerMessage<GadgetMessageIdentifier>* mb =
    new GadgetContainerMessage<GadgetMessageIdentifier>();

  mb->getObjectPtr()->id = GADGET_MESSAGE_ISMRMRD_ACQUISITION;

  mb->cont(m1);
  return controller_->output_ready(mb);
	 */
	throw std::runtime_error("Not implemented.");
	TODO: // Let's try a better implementation.
}

GADGET_FACTORY_DECLARE(AcquisitionFinishGadget)
