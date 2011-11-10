#include "ImageFinishGadget.h"

int ImageFinishGadget
::process(GadgetContainerMessage<GadgetMessageImage>* m1,
	  GadgetContainerMessage< hoNDArray< ACE_UINT16 > >* m2)
{
  if (!controller_) {
    ACE_DEBUG( (LM_DEBUG, 
		ACE_TEXT("Cannot return result to controller, no controller set")) );
    return -1;
  }

  bool is_last_image = (m1->getObjectPtr()->flags & GADGET_FLAG_LAST_IMAGE);

  GadgetContainerMessage<GadgetMessageIdentifier>* mb =
    new GadgetContainerMessage<GadgetMessageIdentifier>();

  mb->getObjectPtr()->id = GADGET_MESSAGE_IMAGE;

  mb->cont(m1);

  int ret =  controller_->output_ready(mb);

  int ret2 = GADGET_OK;
  if (is_last_image) {
	  GadgetContainerMessage<GadgetMessageIdentifier>* mb2 =
	     new GadgetContainerMessage<GadgetMessageIdentifier>();

	  mb2->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;
	  ret2 = controller_->output_ready(mb2);
  }

  if ( (ret < 0) || (ret2 < 0) ) {
	  return GADGET_FAIL;
  }

  return GADGET_OK;
}

GADGET_FACTORY_DECLARE(ImageFinishGadget)
