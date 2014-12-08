#ifndef FFTGADGET_H
#define FFTGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class EXPORTGADGETSMRICORE FFTGadget : 
  public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(FFTGadget)
	
	protected:
      virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1,
			   GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);      
    };
}
#endif //FFTGADGET_H
