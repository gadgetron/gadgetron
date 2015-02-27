#ifndef IMAGEARRAYSPLIT_H
#define IMAGEARRAYSPLIT_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include "mri_core_data.h"

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
