#ifndef FATWATER_H
#define FATWATER_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_fatwater_export.h"

#include "mri_core_data.h"

namespace Gadgetron{

  class EXPORTGADGETFATWATER FatWaterGadget : 
  public Gadget1<IsmrmrdImageArray>
    {
    public:
      GADGET_DECLARE(FatWaterGadget)
      FatWaterGadget();
	
    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdImageArray>* m1);
    };
}
#endif //FATWATER_H
