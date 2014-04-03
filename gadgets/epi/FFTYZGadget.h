#ifndef FFTYZGADGET_H
#define FFTYZGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_epi_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class  EXPORTGADGETS_EPI FFTYZGadget : 
  public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
  {
    public:
      FFTYZGadget();
      virtual ~FFTYZGadget();
	
    protected:
      virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1,
      		           GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);      
  };
}
#endif //FFTYZGADGET_H
