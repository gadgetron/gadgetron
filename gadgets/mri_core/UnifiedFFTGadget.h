#ifndef UNIFIEDFFTGADGET_H
#define UNIFIEDFFTGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include "mri_core_data.h"

namespace Gadgetron{

  class EXPORTGADGETSMRICORE UnifiedFFTGadget : 
  public Gadget1<IsmrmrdReconData>
    {
    public:
      GADGET_DECLARE(UnifiedFFTGadget);

      UnifiedFFTGadget();

    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdReconData>* m1);
      long long image_counter_;
    };

  
}
#endif //UNIFIEDFFTGADGET_H
