/*
 * GrappaUnmixingGadget.h
 *
 *  Created on: Dec 15, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GRAPPAUNMIXINGGADGET_H_
#define GRAPPAUNMIXINGGADGET_H_

#include <Gadget.h>

#include "hoNDArray.h"
#include <complex>
#include "GadgetMRIHeaders.h"
#include "GrappaWeights.h"

struct GrappaUnmixingJob
{
	GrappaWeights<float>* weights_;
};

class GrappaUnmixingGadget: public Gadget3<GrappaUnmixingJob, GadgetMessageImage, hoNDArray<std::complex<float> > > {
public:
	GADGET_DECLARE(GrappaUnmixingGadget);

	GrappaUnmixingGadget();
	virtual ~GrappaUnmixingGadget();
protected:
	virtual int process(GadgetContainerMessage<GrappaUnmixingJob>* m1,
			GadgetContainerMessage<GadgetMessageImage>* m2, GadgetContainerMessage<hoNDArray<std::complex<float> > >* m3);

};

#endif /* GRAPPAUNMIXINGGADGET_H_ */
