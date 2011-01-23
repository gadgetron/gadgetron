#ifndef IMAGEFINISHGADGET_H
#define IMAGEFINISHGADGET_H

#include "Gadget.h"
#include "NDArray.h"
#include "GadgetMRIHeaders.h"
#include "GadgetStreamController.h"

#include <complex>

class ImageFinishGadget : 
public Gadget2<GadgetMessageImage,NDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(ImageFinishGadget);

 protected:
  virtual int process(GadgetContainerMessage<GadgetMessageImage>* m1,
		      GadgetContainerMessage< NDArray< std::complex<float> > >* m2);
  
};


#endif //IMAGEFINISHGADGET_H
