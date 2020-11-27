#include "AcquisitionPassthroughGadget.h"

namespace Gadgetron{
int AcquisitionPassthroughGadget
::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
  //It is enough to put the first one, since they are linked
  if (this->next()->putq(m1) == -1) {
    m1->release();
    GERROR("AcquisitionPassthroughGadget::process, passing data on to next gadget");
    return -1;
  }

  return 0;
}
GADGET_FACTORY_DECLARE(AcquisitionPassthroughGadget)
}


