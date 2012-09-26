#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd.h"
#include "gadgetroncore_export.h"

#include <complex>

class EXPORTGADGETSCORE RemoveROOversamplingGadget :
public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(RemoveROOversamplingGadget);
  
 protected:
  virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
		      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
};

