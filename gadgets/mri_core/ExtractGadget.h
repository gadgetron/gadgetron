#pragma once
#include "Gadget.h"
#include "hoNDArray.h"

#include <bitset>
#include <complex>

namespace Gadgetron {

    class ExtractGadget : public Core::ChannelGadget<mrd::Image<std::complex<float>>>

    {

    public:
        ExtractGadget(const Core::Context& context, const Core::GadgetProperties& props);

    protected:
        NODE_PROPERTY(
            extract_mask, std::bitset<4>, "(DEPRECATED) Extract mask, bitmask MAG=1, REAL=2, IMAG=4, PHASE=8", 0);
        NODE_PROPERTY(extract_magnitude, bool, "Extract absolute value", true);
        NODE_PROPERTY(extract_real, bool, "Extract real components", false);
        NODE_PROPERTY(extract_imag, bool, "Extract imaginary component", false);
        NODE_PROPERTY(extract_phase, bool, "Extract phase", false);
        NODE_PROPERTY(real_imag_offset, float, "Offset to add to real and imag images", 0.0f);

    public:
        void process(Core::InputChannel<mrd::Image<std::complex<float>>>& in, Core::OutputChannel& out) override;

    protected:
        std::set<mrd::ImageType> image_types;
    };
}
