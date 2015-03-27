#ifndef FloatToFixPointAttribGadget_H_
#define FloatToFixPointAttribGadget_H_

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/meta.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>

namespace Gadgetron
{

    /**
    * This Gadget converts float values to fix point integer format.
    *
    * How the conversion is done will depend on the image type:
    * Magnitude images: Values above 4095 will be clamped.
    * Real or Imag: Values below -2048 and above 2047 will be clamped. Zero will be 2048.
    * Phase: -pi will be 0, +pi will be 4095.
    *
    */

    template <typename T> 
    class EXPORTGADGETSMRICORE FloatToFixPointAttribGadget:public Gadget3<ISMRMRD::ImageHeader, hoNDArray< float >, ISMRMRD::MetaContainer >
    {
    public:

        GADGET_DECLARE(FloatToFixPointAttribGadget);

        FloatToFixPointAttribGadget();
        virtual ~FloatToFixPointAttribGadget();

    protected:
        GADGET_PROPERTY(max_intensity, T, "Maximum intensity value", std::numeric_limits<T>::max() );
        GADGET_PROPERTY(min_intensity, T, "Minimal intensity value", std::numeric_limits<T>::min());
        GADGET_PROPERTY(intensity_offset, T, "Intensity offset", 0);

        T max_intensity_value_;
        T min_intensity_value_;
        T intensity_offset_value_;

        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3);
    };

    class EXPORTGADGETSMRICORE FloatToShortAttribGadget :public FloatToFixPointAttribGadget < short > 
    {
    public:
        GADGET_DECLARE(FloatToShortAttribGadget);

        FloatToShortAttribGadget();
        virtual ~FloatToShortAttribGadget();
    };

    class EXPORTGADGETSMRICORE FloatToUShortAttribGadget :public FloatToFixPointAttribGadget < unsigned short >
    {
    public:
        GADGET_DECLARE(FloatToUShortAttribGadget);

        FloatToUShortAttribGadget();
        virtual ~FloatToUShortAttribGadget();
    };

    class EXPORTGADGETSMRICORE FloatToIntAttribGadget :public FloatToFixPointAttribGadget < int >
    {
    public:
        GADGET_DECLARE(FloatToIntAttribGadget);

        FloatToIntAttribGadget();
        virtual ~FloatToIntAttribGadget();
    };

    class EXPORTGADGETSMRICORE FloatToUIntAttribGadget :public FloatToFixPointAttribGadget < unsigned int >
    {
    public:
        GADGET_DECLARE(FloatToUIntAttribGadget);

        FloatToUIntAttribGadget();
        virtual ~FloatToUIntAttribGadget();
    };
}

#endif /* FloatToFixPointAttribGadget_H_ */
