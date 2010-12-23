#include "AcquisitionPassthroughGadget.h"

int AcquisitionPassthroughGadget
::process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
	  GadgetContainerMessage< NDArray< std::complex<float> > >* m2)
{
  ACE_DEBUG( (LM_DEBUG, 
	      ACE_TEXT("AcquisitionPassthroughGadget::process called\n")) );


  //It is enough to put the first one, since they are linked
  if (this->next()->putq(m1) == -1) {
    m1->release();
    ACE_ERROR_RETURN( (LM_ERROR,
		       ACE_TEXT("%p\n"),
		       ACE_TEXT("AcquisitionPassthroughGadget::process, passing data on to next gadget")),
		      -1);
  }

  return 0;
}
