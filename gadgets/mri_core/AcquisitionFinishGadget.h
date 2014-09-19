#ifndef ACQUISITIONFINISHGADGET_H
#define ACQUISITIONFINISHGADGET_H

#include "Gadget.h"
#include "NDArray.h"
#include "GadgetMRIHeaders.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class EXPORTGADGETSMRICORE AcquisitionFinishGadget : 
  public Gadget2<ISMRMRD::AcquisitionHeader, NDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(AcquisitionFinishGadget);
      
    protected:
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< NDArray< std::complex<float> > >* m2);
    };
}

#endif //ACQUISITIONFINISHGADGET_H
