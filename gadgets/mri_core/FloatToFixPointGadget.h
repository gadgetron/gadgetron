#ifndef FloatToFixPointGadget_H_
#define FloatToFixPointGadget_H_

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
    class EXPORTGADGETSMRICORE FloatToFixPointGadget:public Gadget2<ISMRMRD::ImageHeader, hoNDArray< float > >
    {
    public:

        GADGET_DECLARE(FloatToFixPointGadget);

        FloatToFixPointGadget();
        virtual ~FloatToFixPointGadget();

    protected:
        GADGET_PROPERTY(max_intensity, T, "Maximum intensity value", std::numeric_limits<T>::max() );
        GADGET_PROPERTY(min_intensity, T, "Minimal intensity value", std::numeric_limits<T>::min());
        GADGET_PROPERTY(intensity_offset, T, "Intensity offset", 0);

        T max_intensity_value_;
        T min_intensity_value_;
        T intensity_offset_value_;

        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2);
    };

    class EXPORTGADGETSMRICORE FloatToShortGadget :public FloatToFixPointGadget < short > 
    {
    public:
        GADGET_DECLARE(FloatToShortGadget);

        FloatToShortGadget();
        virtual ~FloatToShortGadget();
    };

    class EXPORTGADGETSMRICORE FloatToUShortGadget :public FloatToFixPointGadget < unsigned short >
    {
    public:
        GADGET_DECLARE(FloatToUShortGadget);

        FloatToUShortGadget();
        virtual ~FloatToUShortGadget();
    };

    class EXPORTGADGETSMRICORE FloatToIntGadget :public FloatToFixPointGadget < int >
    {
    public:
        GADGET_DECLARE(FloatToIntGadget);

        FloatToIntGadget();
        virtual ~FloatToIntGadget();
    };

    class EXPORTGADGETSMRICORE FloatToUIntGadget :public FloatToFixPointGadget < unsigned int >
    {
    public:
        GADGET_DECLARE(FloatToUIntGadget);

        FloatToUIntGadget();
        virtual ~FloatToUIntGadget();
    };
}

#endif /* FloatToFixPointGadget_H_ */
