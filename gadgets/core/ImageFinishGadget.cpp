#include "ImageFinishGadget.h"

int ImageFinishGadget
::process(GadgetContainerMessage<GadgetMessageImage>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
  if (!controller_) {
    ACE_DEBUG( (LM_DEBUG, 
		ACE_TEXT("Cannot return result to controller, no controller set")) );
    return -1;
  }

  //Test of permute
  std::vector<unsigned int> dim_order;
  dim_order.push_back(2);
  dim_order.push_back(1);

  if (m2->getObjectPtr()->permute(dim_order) < 0) {
    GADGET_DEBUG1("Permute of array failed\n");
  }

  GadgetContainerMessage<GadgetMessageIdentifier>* mb =
    new GadgetContainerMessage<GadgetMessageIdentifier>();

  mb->getObjectPtr()->id = GADGET_MESSAGE_IMAGE;

  mb->cont(m1);

  return controller_->output_ready(mb);

}

GADGET_FACTORY_DECLARE(ImageFinishGadget)
