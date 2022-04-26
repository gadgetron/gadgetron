#ifndef ACQUISITIONPASSTHROUGHGADGET_H
#define ACQUISITIONPASSTHROUGHGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_debugging_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class EXPORTGADGETSDEBUGGING AcquisitionPassthroughGadget : 
  public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(AcquisitionPassthroughGadget);
      
    protected:
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
    };
}
#endif //ACQUISITIONPASSTHROUGHGADGET_H
