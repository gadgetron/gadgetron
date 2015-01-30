#pragma once
#include "Gadget.h"
#include "hoNDArray.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

#include <boost/make_shared.hpp>

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPUGADGET__)
#define EXPORTGPUGADGET __declspec(dllexport)
#else
#define EXPORTGPUGADGET __declspec(dllimport)
#endif
#else
#define EXPORTGPUGADGET
#endif


namespace Gadgetron{

  class EXPORTGPUGADGET cuFFTGadget :
  public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(cuFFTGadget)

	protected:
      virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1,
			   GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
    };
}
