#ifndef FLOATTOUSHORTGADGET_H_
#define FLOATTOUSHORTGADGET_H_

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>

namespace Gadgetron{
  
  /**
   * This Gadget converts float values to unsigned unsigned short int format.
   *
   * How the conversion is done will depend on the image type:
   * Magnitude images: Values above 4095 will be clamped.
   * Real or Imag: Values below -2048 and above 2047 will be clamped. Zero will be 2048.
   * Phase: -pi will be 0, +pi will be 4095.
   *
   */
  class EXPORTGADGETSMRICORE FloatToUShortGadget:
  public Gadget2<ISMRMRD::ImageHeader,hoNDArray< float > >
    {
    public:

	  FloatToUShortGadget();
      virtual ~FloatToUShortGadget();
      
    protected:
      virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
			  GadgetContainerMessage< hoNDArray< float > >* m2);      
    };
}

#endif /* FLOATTOUSHORTGADGET_H_ */
