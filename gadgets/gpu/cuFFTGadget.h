#pragma once
#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class EXPORTGADGETSMRICORE cuFFTGadget :
  public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(cuFFTGadget)

	protected:
      virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1,
			   GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
    };
}
