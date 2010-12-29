#include "ImageFinishGadget.h"

int ImageFinishGadget
::process(GadgetContainerMessage<GadgetMessageImage>* m1,
	  GadgetContainerMessage< NDArray< std::complex<float> > >* m2)
{
  if (!controller_) {
    ACE_DEBUG( (LM_DEBUG, 
		ACE_TEXT("Cannot return result to controller, no controller set")) );
    return -1;
  }
  
  GadgetContainerMessage<GadgetMessageIdentifier>* mb =
    new GadgetContainerMessage<GadgetMessageIdentifier>();

  mb->getObjectPtr()->id = GADGET_MESSAGE_IMAGE;

  mb->cont(m1);

  return controller_->output_ready(mb);

}
