/*
 * AutoScaleGadget.h
 *
 *  Created on: Dec 19, 2011
 *      Author: Michael S. Hansen
 */

#ifndef AUTOSCALEGADGET_H_
#define AUTOSCALEGADGET_H_

#include <Gadget.h>
#include "ismrmrd.h"
#include "hoNDArray.h"

class AutoScaleGadget:
public Gadget2<ISMRMRD::ImageHeader,hoNDArray< float > >
{
public:
	GADGET_DECLARE(AutoScaleGadget);
	AutoScaleGadget();
	virtual ~AutoScaleGadget();
protected:
	virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
			GadgetContainerMessage< hoNDArray< float > >* m2);

	unsigned int histogram_bins_;
	std::vector<unsigned int> histogram_;
	float current_scale_;
	float max_value_;

};

#endif /* AUTOSCALEGADGET_H_ */
