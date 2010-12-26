#ifndef ACQUISITIONFINISHGADGET_H
#define ACQUISITIONFINISHGADGET_H

#include "Gadget.h"
#include "NDArray.h"
#include "gadgetheaders.h"
#include "GadgetStreamController.h"

#include <complex>

class AcquisitionFinishGadget : 
public Gadget2<GadgetMessageAcquisition,NDArray< std::complex<float> > >
{
 public:
  AcquisitionFinishGadget(GadgetStreamController* controller) 
    : controller_(controller)
    { }
  virtual int process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
		      GadgetContainerMessage< NDArray< std::complex<float> > >* m2);

 protected:
  GadgetStreamController* controller_;
  
};


#endif //ACQUISITIONFINISHGADGET_H
