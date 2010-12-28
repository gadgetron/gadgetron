#include "FFTGadget.h"
#include "FFT.h"

int FFTGadget::process( GadgetContainerMessage< GadgetMessageImage>* m1,
			GadgetContainerMessage< NDArray< std::complex<float> > >* m2)
{
  FFT<float> ft;

  ft.ifft(m2->getObjectPtr(),0);
  ft.ifft(m2->getObjectPtr(),1);
  ft.ifft(m2->getObjectPtr(),2);

  return this->next()->putq(m1);
}
