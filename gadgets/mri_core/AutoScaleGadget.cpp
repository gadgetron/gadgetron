/*
 * AutoScaleGadget.cpp
 *
 *  Created on: Dec 19, 2011
 *      Author: Michael S. Hansen
 */

#include "AutoScaleGadget.h"

namespace Gadgetron{

AutoScaleGadget::AutoScaleGadget()
	: histogram_bins_(100)
	, current_scale_(1.0)
	, max_value_(2048)
{
}

AutoScaleGadget::~AutoScaleGadget() {
	// TODO Auto-generated destructor stub
}

int AutoScaleGadget::process_config(ACE_Message_Block* mb) {
        max_value_ = max_value.value();
	return GADGET_OK;
}


int AutoScaleGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2)
{
	if (m1->getObjectPtr()->image_type == ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE) { //Only scale magnitude images for now
		float max = 0.0f;
		float* d = m2->getObjectPtr()->get_data_ptr();
		for (unsigned long int i = 0; i < m2->getObjectPtr()->get_number_of_elements(); i++) {
			if (d[i] > max) max = d[i];
		}

		if (histogram_.size() != histogram_bins_) {
			histogram_ = std::vector<size_t>(histogram_bins_);
		}

		for (size_t i = 0; i < histogram_bins_; i++) {
			histogram_[i] = 0;
		}

		for (unsigned long int i = 0; i < m2->getObjectPtr()->get_number_of_elements(); i++) {
			size_t bin = static_cast<size_t>(std::floor((d[i]/max)*histogram_bins_));
			if (bin >= histogram_bins_) {
				bin = histogram_bins_-1;
			}
			histogram_[bin]++;
		}

		//Find 99th percentile
		long long cumsum = 0;
		size_t counter = 0;
		while (cumsum < (0.99*m2->getObjectPtr()->get_number_of_elements())) {
			cumsum += (long long)(histogram_[counter++]);
		}
		max = (counter+1)*(max/histogram_bins_);
		GDEBUG("Max: %f\n",max);

		current_scale_ = max_value_/max;

		for (unsigned long int i = 0; i < m2->getObjectPtr()->get_number_of_elements(); i++) {
			d[i] *= current_scale_;
		}
	}

	if (this->next()->putq(m1) < 0) {
		GDEBUG("Failed to pass on data to next Gadget\n");
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(AutoScaleGadget)

}
