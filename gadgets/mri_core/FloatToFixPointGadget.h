#pragma once

#include "Node.h"

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

    template <typename T, typename Base >
    class FloatToFixPointGadget: public Core::ChannelGadget<mrd::Image<float>>
    {
    public:

        using Core::ChannelGadget<mrd::Image<float>>::ChannelGadget;

        ~FloatToFixPointGadget() override = default ;

        void process(Core::InputChannel<mrd::Image<float>>& input, Core::OutputChannel& output) override;
    };

    class FloatToShortGadget :public FloatToFixPointGadget < short,FloatToShortGadget >
    {
    public:
        using FloatToFixPointGadget<short,FloatToShortGadget>::FloatToFixPointGadget;
        NODE_PROPERTY(max_intensity, short, "Maximum intensity value", std::numeric_limits<short >::max() );
        NODE_PROPERTY(min_intensity, short , "Minimal intensity value", std::numeric_limits<short >::min());
        NODE_PROPERTY(intensity_offset, short , "Intensity offset", 0);
        ~FloatToShortGadget() override = default;
    };
    class FloatToUShortGadget :public FloatToFixPointGadget < unsigned short,FloatToUShortGadget >
    {
    public:
        using FloatToFixPointGadget<unsigned short,FloatToUShortGadget>::FloatToFixPointGadget;
        NODE_PROPERTY(max_intensity, short, "Maximum intensity value", 4095 );
        NODE_PROPERTY(min_intensity, short , "Minimal intensity value", 0);
        NODE_PROPERTY(intensity_offset, short , "Intensity offset", 2048);
        ~FloatToUShortGadget() override = default;
    };
    class FloatToUIntGadget :public FloatToFixPointGadget < unsigned int,FloatToUIntGadget >
    {
    public:
        using FloatToFixPointGadget<unsigned int,FloatToUIntGadget>::FloatToFixPointGadget;
        NODE_PROPERTY(max_intensity, int, "Maximum intensity value", 4095 );
        NODE_PROPERTY(min_intensity, int , "Minimal intensity value", 0);
        NODE_PROPERTY(intensity_offset, int , "Intensity offset", 2048);
        ~FloatToUIntGadget() override = default;
    };
    class FloatToIntGadget :public FloatToFixPointGadget < int,FloatToIntGadget >
    {
    public:
        using FloatToFixPointGadget<int,FloatToIntGadget>::FloatToFixPointGadget;
        NODE_PROPERTY(max_intensity, int, "Maximum intensity value", 4095 );
        NODE_PROPERTY(min_intensity, int , "Minimal intensity value", 0);
        NODE_PROPERTY(intensity_offset, int , "Intensity offset", 2048);
        ~FloatToIntGadget() override = default;
    };





}

