#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"

#include <complex>

class EXPORTGADGETSGRAPPA RemoveROOversamplingGadget : 
public Gadget2<GadgetMessageAcquisition,hoNDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(RemoveROOversamplingGadget);
  
 protected:
  virtual int process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
		      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
};

