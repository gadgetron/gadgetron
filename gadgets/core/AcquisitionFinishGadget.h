#ifndef ACQUISITIONFINISHGADGET_H
#define ACQUISITIONFINISHGADGET_H

#include "gadgetron_core_export.h"
#include "Gadget.h"
#include "NDArray.h"
#include "ismrmrd.h"
#include "GadgetMRIHeaders.h"

#include <complex>

namespace Gadgetron{

class EXPORTGADGETSCORE AcquisitionFinishGadget : 
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
