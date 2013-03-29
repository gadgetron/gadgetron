#include "GadgetIsmrmrdReadWrite.h"
#include "FFTGadget.h"
#include "hoFFT.h"

namespace Gadgetron{

  int FFTGadget::process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    hoFFT<float>::instance()->ifft(m2->getObjectPtr(),0);
    hoFFT<float>::instance()->ifft(m2->getObjectPtr(),1);
    hoFFT<float>::instance()->ifft(m2->getObjectPtr(),2);
    
    if (this->next()->putq(m1) < 0) {
      return GADGET_FAIL;
    }
    
    return GADGET_OK;    
  }
  
  GADGET_FACTORY_DECLARE(FFTGadget)
}
