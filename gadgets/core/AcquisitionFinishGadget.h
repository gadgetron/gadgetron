#ifndef ACQUISITIONFINISHGADGET_H
#define ACQUISITIONFINISHGADGET_H

#include "gadgetron_export.h"
#include "Gadget.h"
#include "NDArray.h"
#include "GadgetMRIHeaders.h"

#include <complex>

class EXPORTGADGETSCORE AcquisitionFinishGadget : 
public Gadget2<GadgetMessageAcquisition,NDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(AcquisitionFinishGadget);
  
 protected:
  virtual int process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
		      GadgetContainerMessage< NDArray< std::complex<float> > >* m2);

};


#endif //ACQUISITIONFINISHGADGET_H
