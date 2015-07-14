#ifndef IMAGESPLITSLICES_H
#define IMAGESPLITSLICES_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include "mri_core_data.h"

namespace Gadgetron{

  class EXPORTGADGETSMRICORE ImageSplitSlicesGadget :
  public Gadget2<ISMRMRD::ImageHeader,hoNDArray< std::complex<float> > >
  {
    public:
      GADGET_DECLARE(ImageSplitSlicesGadget)
      ImageSplitSlicesGadget();
	
    protected:
      virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
    };
}
#endif //IMAGESPLITSLICES_H
