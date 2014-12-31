#include "GadgetIsmrmrdReadWrite.h"
#include "ImageFinishGadget.h"
namespace Gadgetron{
template <typename T>
int ImageFinishGadget<T>
::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
	  GadgetContainerMessage< hoNDArray< T > >* m2)
{
  if (!this->controller_) {
    ACE_DEBUG( (LM_DEBUG, 
		ACE_TEXT("Cannot return result to controller, no controller set")) );
    return -1;
  }

  GadgetContainerMessage<GadgetMessageIdentifier>* mb =
    new GadgetContainerMessage<GadgetMessageIdentifier>();

  switch (sizeof(T)) {
  case 2: //Unsigned short
	  mb->getObjectPtr()->id = GADGET_MESSAGE_IMAGE_REAL_USHORT;
	  break;
  case 4: //Float
	  mb->getObjectPtr()->id = GADGET_MESSAGE_IMAGE_REAL_FLOAT;
	  break;
  case 8: //Complex float
	  mb->getObjectPtr()->id = GADGET_MESSAGE_IMAGE_CPLX_FLOAT;
	  break;
  default:
	  GDEBUG("Wrong data size detected: %d\n", sizeof(T));
	  mb->release();
	  m1->release();
	  return GADGET_FAIL;
  }

  mb->cont(m1);

  int ret =  this->controller_->output_ready(mb);

  if ( (ret < 0) ) {
	  GDEBUG("Failed to return massage to controller\n");
	  return GADGET_FAIL;
  }

  return GADGET_OK;
}

//Declare factories for the various template instances
GADGET_FACTORY_DECLARE(ImageFinishGadgetFLOAT);
GADGET_FACTORY_DECLARE(ImageFinishGadgetUSHORT);
GADGET_FACTORY_DECLARE(ImageFinishGadgetCPLX);
}
