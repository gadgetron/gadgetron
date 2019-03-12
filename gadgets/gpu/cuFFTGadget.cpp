
#include "cuFFTGadget.h"
#include "cuNDFFT.h"

namespace Gadgetron{

  int cuFFTGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  	{
	  hoNDArray<complext<float>> * tmp = (hoNDArray<complext<float>>*) m2->getObjectPtr();
	  cuNDArray< complext<float> > cu_data(*tmp);

	  cu_data.squeeze();
	  cuNDFFT<float>::instance()->ifft(&cu_data);
	  cu_data.to_host(tmp);

    if (this->next()->putq(m1) < 0) {
      return GADGET_FAIL;
    }

    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(cuFFTGadget)
}
