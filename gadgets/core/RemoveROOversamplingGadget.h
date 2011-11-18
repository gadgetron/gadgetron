#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "gadgetroncore_export.h"

#include <complex>

class EXPORTGADGETSCORE RemoveROOversamplingGadget :
public Gadget2<GadgetMessageAcquisition,hoNDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(RemoveROOversamplingGadget);
  
 protected:
  virtual int process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
		      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
};

