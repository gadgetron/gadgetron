#ifndef IMAGEFINISHGADGET_H
#define IMAGEFINISHGADGET_H

#include "Gadget.h"
#include "NDArray.h"
#include "gadgetheaders.h"
#include "GadgetStreamController.h"

#include <complex>

class ImageFinishGadget : 
public Gadget2<GadgetMessageImage,NDArray< std::complex<float> > >
{
 public:
  ImageFinishGadget(GadgetStreamController* controller) 
    : controller_(controller)
    { }
  virtual int process(GadgetContainerMessage<GadgetMessageImage>* m1,
		      GadgetContainerMessage< NDArray< std::complex<float> > >* m2);

 protected:
  GadgetStreamController* controller_;
  
};


#endif //IMAGEFINISHGADGET_H
