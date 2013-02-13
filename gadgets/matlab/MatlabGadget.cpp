#include "MatlabGadget.h"


int AcquisitionMatlabGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

	//GADGET_DEBUG1("Data passing MATLAB Gadget\n");

	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}



	return GADGET_OK;
}

int ImageMatlabGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}


	return GADGET_OK;
}


GADGET_FACTORY_DECLARE(AcquisitionMatlabGadget)
GADGET_FACTORY_DECLARE(ImageMatlabGadget)
