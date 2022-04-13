/**
    \brief  Performs denoising with non-local means and non-local Bayes
    \author Original: Thomas Sangild Sorensen
    \author PureGadget Conversion: David Christoffer Hansen
    \test   Tested by:epi_2d.cfg and others
*/

#pragma once
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "gadgetron_mricore_export.h"
#include "hoNDArray.h"

#include <bitset>
#include <complex>
#include <ismrmrd/ismrmrd.h>

namespace Gadgetron {

    class ExtractGadget : public Core::ChannelGadget<Core::Image<std::complex<float>>>

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
        void process(Core::InputChannel<Core::Image<std::complex<float>>>& in, Core::OutputChannel& out) override;

    protected:
        std::set<ISMRMRD::ISMRMRD_ImageTypes> image_types;
    };
}
