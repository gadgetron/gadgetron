#ifndef FFTGADGET_H
#define FFTGADGET_H

#include "gadgetron_core_export.h"
#include "Gadget.h"
#include "ismrmrd.h"
#include "hoNDArray.h"

#include <complex>

namespace Gadgetron{

class EXPORTGADGETSCORE FFTGadget : 
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
