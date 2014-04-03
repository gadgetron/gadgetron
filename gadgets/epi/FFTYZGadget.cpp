#include "GadgetIsmrmrdReadWrite.h"
#include "FFTYZGadget.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"

namespace Gadgetron{

  FFTYZGadget::FFTYZGadget() {}
  FFTYZGadget::~FFTYZGadget() {}

  int FFTYZGadget::process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1,
			    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {

    // FFT along 1st two dimensions (y and z)
    hoNDFFT<float>::instance()->ifft(m2->getObjectPtr(),1);
    hoNDFFT<float>::instance()->ifft(m2->getObjectPtr(),2);

    if (this->next()->putq(m1) < 0) {
      return GADGET_FAIL;
    }
    
    return GADGET_OK;    
  }
  
  GADGET_FACTORY_DECLARE(FFTYZGadget)
}
