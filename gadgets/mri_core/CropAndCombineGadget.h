#ifndef CROPANDCOMBINEGADGET_H
#define CROPANDCOMBINEGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class CropAndCombineGadget :
  public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    protected:
      virtual int process( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
			   GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
    };
}

#endif //CROPANDCOMBINEGADGET_H
