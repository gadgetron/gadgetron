#ifndef CROPANDCOMBINEGADGET_H
#define CROPANDCOMBINEGADGET_H

#include "gadgetron_core_export.h"
#include "Gadget.h"
#include "ismrmrd.h"
#include "hoNDArray.h"

#include <complex>

namespace Gadgetron{

class EXPORTGADGETSCORE CropAndCombineGadget : 
public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(CropAndCombineGadget);

 protected:
  virtual int process( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
		       GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
		     
};
}
#endif //CROPANDCOMBINEGADGET_H
