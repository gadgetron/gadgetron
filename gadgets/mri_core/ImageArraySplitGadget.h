#ifndef IMAGEARRAYSPLIT_H
#define IMAGEARRAYSPLIT_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include "mri_core_data.h"

// this gadget will not copy the waveform data down the chain for every image

namespace Gadgetron{

  class EXPORTGADGETSMRICORE ImageArraySplitGadget : 
  public Gadget1<IsmrmrdImageArray>
    {
    public:
      GADGET_DECLARE(ImageArraySplitGadget)
      ImageArraySplitGadget();
	
    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdImageArray>* m1);
    };
}
#endif //IMAGEARRAYSPLIT_H
