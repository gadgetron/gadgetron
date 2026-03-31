/**
    \brief  Apply extra scaling to the readout, to compensate for incorrect noise scaling.
    \author Original: Hui Xue
*/

#pragma once

#include "Node.h"
#include "Types.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "Types.h"
#include "hoNDArray.h"

namespace Gadgetron {
class ReadoutScalingGadget : public Core::ChannelGadget<Core::Acquisition>
{
    public:
        using Core::ChannelGadget<Core::Acquisition>::ChannelGadget;

        ReadoutScalingGadget(const Core::Context& context, const Core::GadgetProperties& props);
        ~ReadoutScalingGadget() override = default;

        void process(Core::InputChannel<Core::Acquisition>& input, Core::OutputChannel& output) override;

    protected:
        NODE_PROPERTY(noise_scaling_factor, double, "For 256 noise lines, 512 sample each line", 17179869184);
};

} // namespace Gadgetron