#ifndef COMBINEGADGET_H
#define COMBINEGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "hoArmadillo.h"
#include "gadgetron_epi_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron{
  
  class  EXPORTGADGETS_EPI CombineGadget : 
  public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(CombineGadget);
      
    protected:
      virtual int process( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
			   GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);     
    };
}

#endif //COMBINEGADGET_H
