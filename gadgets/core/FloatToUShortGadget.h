/*
 * FloatToUShortGadget.h
 *
 *  Created on: Nov 26, 2011
 *      Author: hansenms
 */

#ifndef FLOATTOUSHORTGADGET_H_
#define FLOATTOUSHORTGADGET_H_

#include <Gadget.h>
#include <hoNDArray.h>
#include "GadgetMRIHeaders.h"
#include "gadgetroncore_export.h"

/**
 * This Gadget converts float values to unsigned unsigned short int format.
 *
 * How the conversion is done will depend on the image type:
 * Magnitude images: Values above 4095 will be clamped.
 * Real or Imag: Values below -2048 and above 2047 will be clamped. Zero will be 2048.
 * Phase: -pi will be 0, +pi will be 4095.
 *
 */
class EXPORTGADGETSCORE FloatToUShortGadget:
public Gadget2<GadgetMessageImage,hoNDArray< float > >
{
public:
	GADGET_DECLARE(FloatToUShortGadget);
	FloatToUShortGadget();
	virtual ~FloatToUShortGadget();

protected:
	virtual int process(GadgetContainerMessage<GadgetMessageImage>* m1,
			GadgetContainerMessage< hoNDArray< float > >* m2);

};




#endif /* FLOATTOUSHORTGADGET_H_ */
