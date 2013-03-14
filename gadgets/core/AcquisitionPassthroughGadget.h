#ifndef ACQUISITIONPASSTHROUGHGADGET_H
#define ACQUISITIONPASSTHROUGHGADGET_H

#include "gadgetroncore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd.h"

#include <complex>
namespace Gadgetron{
class EXPORTGADGETSCORE AcquisitionPassthroughGadget : 
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
