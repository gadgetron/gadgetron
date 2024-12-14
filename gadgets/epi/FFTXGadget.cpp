#include "FFTXGadget.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"

namespace Gadgetron{

  FFTXGadget::FFTXGadget() {}
  FFTXGadget::~FFTXGadget() {}

  int FFTXGadget::process(GadgetContainerMessage< mrd::Acquisition>* m1)
  {

    // FFT along 1st dimensions (x)
    if(buf_.get_size(0)!= m1->getObjectPtr()->data.get_size(0))
    {
        buf_ = m1->getObjectPtr()->data;
    }

    hoNDFFT<float>::instance()->fft1c( m1->getObjectPtr()->data, r_, buf_);
    memcpy(m1->getObjectPtr()->data.begin(), r_.begin(), r_.get_number_of_bytes());

    if (this->next()->putq(m1) < 0)
    {
      return GADGET_FAIL;
    }
    
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(FFTXGadget)
}
