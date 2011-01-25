#ifndef CROPANDCOMBINEGADGET_H
#define CROPANDCOMBINEGADGET_H

#include <complex>

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "NDArray.h"

class CropAndCombineGadget : 
public Gadget2<GadgetMessageImage, NDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(CropAndCombineGadget);

 protected:
  virtual int process( GadgetContainerMessage<GadgetMessageImage>* m1,
		       GadgetContainerMessage< NDArray< std::complex<float> > >* m2);
		     
};

#endif //CROPANDCOMBINEGADGET_H
