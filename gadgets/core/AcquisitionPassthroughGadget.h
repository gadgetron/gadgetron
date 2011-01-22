#ifndef ACQUISITIONPASSTHROUGHGADGET_H
#define ACQUISITIONPASSTHROUGHGADGET_H

#include "Gadget.h"
#include "NDArray.h"
#include "gadgetheaders.h"

#include <complex>

class AcquisitionPassthroughGadget : 
public Gadget2<GadgetMessageAcquisition,NDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(AcquisitionPassthroughGadget);
  
 protected:
  virtual int process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
		      GadgetContainerMessage< NDArray< std::complex<float> > >* m2);
};


#endif //ACQUISITIONPASSTHROUGHGADGET_H
