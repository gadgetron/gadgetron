/*
 * ExtractMagnitudeGadget.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: Michael S. Hansen
 */

#include "ExtractMagnitudeGadget.h"


int ExtractMagnitudeGadget::process(GadgetContainerMessage<GadgetMessageImage> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{

	GadgetContainerMessage<hoNDArray< ACE_UINT16 > > *m3 =
			new GadgetContainerMessage<hoNDArray< ACE_UINT16 > >();

	boost::shared_ptr< std::vector<unsigned int> > dims = m2->getObjectPtr()->get_dimensions();

	if (!m3->getObjectPtr()->create(dims.get())) {
		GADGET_DEBUG1("Unable to create unsigned short storage in Extract Magnitude Gadget");
		return GADGET_FAIL;
	}

	std::complex<float>* src = m2->getObjectPtr()->get_data_ptr();
	ACE_UINT16* dst = m3->getObjectPtr()->get_data_ptr();

	for (unsigned long i = 0; i < m3->getObjectPtr()->get_number_of_elements(); i++) {
		float pix_val = abs(src[i]);
		if (pix_val > 4095) pix_val = 4095;
		dst[i] = static_cast<ACE_UINT16>(pix_val);
	}

	m1->cont(m3);
	m2->release();

	m1->getObjectPtr()->image_format = GADGET_IMAGE_REAL_UNSIGNED_SHORT;
	m1->getObjectPtr()->image_type = GADGET_IMAGE_MAGNITUDE;

	if (this->next()->putq(m1) == -1) {
		m1->release();
		GADGET_DEBUG1("Unable to put unsigned short magnitude image on next gadgets queue");
		return GADGET_FAIL;
	}



	return GADGET_OK;
}


GADGET_FACTORY_DECLARE(ExtractMagnitudeGadget)

