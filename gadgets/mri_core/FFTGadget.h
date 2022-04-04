#pragma once
#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include "mri_core_data.h"

namespace Gadgetron{

  class EXPORTGADGETSMRICORE FFTGadget : 
  public Gadget1<IsmrmrdReconData>
    {
    public:
      GADGET_DECLARE(FFTGadget)
      FFTGadget();
	
    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdReconData>* m1);
      long long image_counter_;      
    };
}
