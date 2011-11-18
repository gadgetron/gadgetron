#ifndef IMAGEFINISHGADGET_H
#define IMAGEFINISHGADGET_H

#include "gadgetroncore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "GadgetStreamController.h"

#include <complex>

class EXPORTGADGETSCORE ImageFinishGadget : 
public Gadget2<GadgetMessageImage,hoNDArray< ACE_UINT16 > >
{
 public:
  GADGET_DECLARE(ImageFinishGadget);

 protected:
  virtual int process(GadgetContainerMessage<GadgetMessageImage>* m1,
		      GadgetContainerMessage< hoNDArray< ACE_UINT16 > >* m2);
  
};


#endif //IMAGEFINISHGADGET_H
