/*
 * AutoScaleGadget.cpp
 *
 *  Created on: Dec 19, 2011
 *      Author: Michael S. Hansen
 */

#include "GadgetIsmrmrdReadWrite.h"
#include "AutoScaleGadget.h"

AutoScaleGadget::AutoScaleGadget()
	: histogram_bins_(100)
	, current_scale_(1.0)
	, max_value_(2048)
{

}

AutoScaleGadget::~AutoScaleGadget() {
	// TODO Auto-generated destructor stub
}


int AutoScaleGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2)
{
	if (m1->getObjectPtr()->image_type == ISMRMRD::TYPE_MAGNITUDE) { //Only scale magnitude images for now
		float max = 0.0f;
		float* d = m2->getObjectPtr()->get_data_ptr();
		for (unsigned long int i = 0; i < m2->getObjectPtr()->get_number_of_elements(); i++) {
			if (d[i] > max) max = d[i];
		}

		if (histogram_.size() != histogram_bins_) {
			histogram_ = std::vector<unsigned int>(histogram_bins_);
		}

		for (unsigned int i = 0; i < histogram_bins_; i++) {
			histogram_[i] = 0;
		}

		for (unsigned long int i = 0; i < m2->getObjectPtr()->get_number_of_elements(); i++) {
			unsigned int bin = static_cast<unsigned int>(floor((d[i]/max)*histogram_bins_));
			if (bin >= histogram_bins_) {
				bin = histogram_bins_-1;
			}
			histogram_[bin]++;
		}

		//Find 99th percentile
		long cumsum = 0;
		unsigned int counter = 0;
		while (cumsum < (0.99*m2->getObjectPtr()->get_number_of_elements())) {
			cumsum += histogram_[counter++];
		}
		max = (counter+1)*(max/histogram_bins_);

		current_scale_ = max_value_/max;

		for (unsigned long int i = 0; i < m2->getObjectPtr()->get_number_of_elements(); i++) {
			d[i] *= current_scale_;
		}
	}

	if (this->next()->putq(m1) < 0) {
		GADGET_DEBUG1("Failed to pass on data to next Gadget\n");
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(AutoScaleGadget)


