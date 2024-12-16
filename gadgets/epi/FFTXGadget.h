#ifndef FFTXGADGET_H
#define FFTXGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"

#include <complex>

namespace Gadgetron{

  class FFTXGadget : 
  public Gadget1<mrd::Acquisition>
  {
    public:
      FFTXGadget();
      virtual ~FFTXGadget();

    protected:
      virtual int process( GadgetContainerMessage< mrd::Acquisition>* m1);

      hoNDArray< std::complex<float> > r_;
      hoNDArray< std::complex<float> > buf_;
  };
}
#endif //FFTXGADGET_H
