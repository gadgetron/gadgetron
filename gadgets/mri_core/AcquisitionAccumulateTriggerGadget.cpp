#include "GadgetIsmrmrdReadWrite.h"
#include "AcquisitionAccumulateTriggerGadget.h"
#include "Gadgetron.h"
#include "mri_core_data.h"

namespace Gadgetron{

  int AcquisitionAccumulateTriggerGadget
  ::process_config(ACE_Message_Block* mb)
  {

    std::string trigger_dimension = *this->get_string_value("trigger_dimension");
  
    GADGET_DEBUG2("TRIGGER DIMENSION IS: %s\n", trigger_dimension.c_str());

    return GADGET_OK;
  }

  int AcquisitionAccumulateTriggerGadget
  ::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
	    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    //It is enough to put the first one, since they are linked
    if (this->next()->putq(m1) == -1) {
      m1->release();
      ACE_ERROR_RETURN( (LM_ERROR,
			 ACE_TEXT("%p\n"),
			 ACE_TEXT("AcquisitionAccumulateTriggerGadget::process, passing data on to next gadget")),
			-1);
    }
    
    return GADGET_OK;
  }

GADGET_FACTORY_DECLARE(AcquisitionAccumulateTriggerGadget)

}


