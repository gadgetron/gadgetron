#ifndef FFTXGADGET_H
#define FFTXGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_epi_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class   EXPORTGADGETS_EPI FFTXGadget : 
  public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
  {
    public:
      FFTXGadget();
      virtual ~FFTXGadget();

    protected:
      virtual int process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader>* m1,
                       GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
  };
}
#endif //FFTXGADGET_H
