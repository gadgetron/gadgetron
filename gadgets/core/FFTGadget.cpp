#include "FFTGadget.h"
#include "FFT.h"

int FFTGadget::process( GadgetContainerMessage< GadgetMessageImage>* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
  FFT<float>::instance()->ifft(m2->getObjectPtr(),0);
  FFT<float>::instance()->ifft(m2->getObjectPtr(),1);
  FFT<float>::instance()->ifft(m2->getObjectPtr(),2);

  if (this->next()->putq(m1) < 0) {
     return GADGET_FAIL;
  }

  return GADGET_OK;

}

GADGET_FACTORY_DECLARE(FFTGadget)
