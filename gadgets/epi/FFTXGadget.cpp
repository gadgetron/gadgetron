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
    if(buf_.get_size(0)!= m2->getObjectPtr()->get_size(0))
    {
        buf_ = *m2->getObjectPtr();
    }

    hoNDFFT<float>::instance()->fft1c( *m2->getObjectPtr(), r_, buf_);
    memcpy(m2->getObjectPtr()->begin(), r_.begin(), r_.get_number_of_bytes());

    if (this->next()->putq(m1) < 0)
    {
      return GADGET_FAIL;
    }
    
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(FFTXGadget)
}
