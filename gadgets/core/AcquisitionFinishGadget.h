#ifndef ACQUISITIONFINISHGADGET_H
#define ACQUISITIONFINISHGADGET_H

#include "gadgetroncore_export.h"
#include "Gadget.h"
#include "NDArray.h"
#include "ismrmrd.h"
#include "GadgetMRIHeaders.h"

#include <complex>

class EXPORTGADGETSCORE AcquisitionFinishGadget : 
public Gadget2<ISMRMRD::AcquisitionHeader, GADGETRON::NDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(AcquisitionFinishGadget);
  
 protected:
  virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
		      GadgetContainerMessage< GADGETRON::NDArray< std::complex<float> > >* m2);

};


#endif //ACQUISITIONFINISHGADGET_H
