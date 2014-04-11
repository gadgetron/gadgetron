#include "GadgetIsmrmrdReadWrite.h"
#include "FFTXGadget.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"

namespace Gadgetron{

  FFTXGadget::FFTXGadget() {}
  FFTXGadget::~FFTXGadget() {}

  int FFTXGadget::process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
			    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {

    // FFT along 1st dimensions (x)
    hoNDFFT<float>::instance()->fft(m2->getObjectPtr(),0);

    if (this->next()->putq(m1) < 0) {
      return GADGET_FAIL;
    }
    
    return GADGET_OK;    
  }
  
  GADGET_FACTORY_DECLARE(FFTXGadget)
}
