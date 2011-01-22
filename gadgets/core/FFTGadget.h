#ifndef FFTGADGET_H
#define FFTGADGET_H

#include "Gadget.h"
#include "gadgetheaders.h"
#include "NDArray.h"
#include <complex>

class FFTGadget : 
public Gadget2<GadgetMessageImage, NDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(FFTGadget)

 protected:
  virtual int process( GadgetContainerMessage< GadgetMessageImage>* m1,
		       GadgetContainerMessage< NDArray< std::complex<float> > >* m2);

};

#endif //FFTGADGET_H
