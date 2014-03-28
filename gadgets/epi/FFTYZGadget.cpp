#include "GadgetIsmrmrdReadWrite.h"
#include "FFTYZGadget.h"
#include "hoNDFFT.h"

namespace Gadgetron{

  int FFTYZGadget::process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1,
			    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    // FFT along ky and kz
    hoNDFFT<float>::instance()->ifft(m2->getObjectPtr(),1);
    hoNDFFT<float>::instance()->ifft(m2->getObjectPtr(),2);
    
    if (this->next()->putq(m1) < 0) {
      return GADGET_FAIL;
    }
    
    return GADGET_OK;    
  }
  
  GADGET_FACTORY_DECLARE(FFTGadget)
}
