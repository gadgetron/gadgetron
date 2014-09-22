#ifndef FloatToUShortAttribGadget_H_
#define FloatToUShortAttribGadget_H_

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/meta.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>

namespace Gadgetron
{

    /**
    * This Gadget converts float values to unsigned unsigned short int format.
    *
    * How the conversion is done will depend on the image type:
    * Magnitude images: Values above 4095 will be clamped.
    * Real or Imag: Values below -2048 and above 2047 will be clamped. Zero will be 2048.
    * Phase: -pi will be 0, +pi will be 4095.
    *
    */

    class EXPORTGADGETSMRICORE FloatToUShortAttribGadget:public Gadget3<ISMRMRD::ImageHeader, hoNDArray< float >, ISMRMRD::MetaContainer >
    {
    public:

        GADGET_DECLARE(FloatToUShortAttribGadget);

        FloatToUShortAttribGadget();
        virtual ~FloatToUShortAttribGadget();

    protected:

        ACE_UINT16 max_intensity_value_;
        ACE_UINT16 intensity_offset_value_;

        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3);
    };
}

#endif /* FloatToUShortAttribGadget_H_ */
