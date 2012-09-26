#ifndef CROPANDCOMBINEGADGET_H
#define CROPANDCOMBINEGADGET_H

#include <complex>

#include "gadgetroncore_export.h"
#include "Gadget.h"
#include "ismrmrd.h"
#include "hoNDArray.h"

class EXPORTGADGETSCORE CropAndCombineGadget : 
public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(CropAndCombineGadget);

 protected:
  virtual int process( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
		       GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
		     
};

#endif //CROPANDCOMBINEGADGET_H
