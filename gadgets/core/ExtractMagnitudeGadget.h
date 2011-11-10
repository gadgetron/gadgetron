/*
 * ExtractMagnitudeGadget.h
 *
 *  Created on: Nov 8, 2011
 *      Author: Michael S. Hansen
 */

#ifndef EXTRACTMAGNITUDEGADGET_H_
#define EXTRACTMAGNITUDEGADGET_H_

#include <Gadget.h>
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "gadgetron_export.h"

#include <complex>

#define MAX_UNSIGNED_SHORT_IMAGE_VALUE

/**
 * This Gadget extracts the magnitude of complex images and converts to an unsigned short int format.
 * Values will be scaled in the range of [0 ... 4095]. Values above 4095 will be clamped at 4095.
 *
 */
class EXPORTGADGETSCORE ExtractMagnitudeGadget:
public Gadget2<GadgetMessageImage,hoNDArray< std::complex<float> > >
{

public:
  GADGET_DECLARE(ExtractMagnitudeGadget);

protected:
	virtual int process(GadgetContainerMessage<GadgetMessageImage>* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
};

#endif /* EXTRACTMAGNITUDEGADGET_H_ */
