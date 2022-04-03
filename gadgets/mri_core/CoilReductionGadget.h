/**
    \brief  
    \author Original: Michael S. Hansen
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Tested by: simple_gre_3d.cfg, gpu_spiral_realtime_deblurring.cfg, gpu_fixed_radial_mode1_realtime.cfg, and others
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"

namespace Gadgetron{
  class CoilReductionGadget : public Core::ChannelGadget<Core::Acquisition> 
    {
      public:
        using Core::ChannelGadget<Core::Acquisition>::ChannelGadget;
        CoilReductionGadget(const Core::Context& context, const Core::GadgetProperties& props);
        ~CoilReductionGadget() override = default;
        void process(Core::InputChannel<Core::Acquisition>& input, Core::OutputChannel& output) override;
      protected:
        NODE_PROPERTY(coil_mask, std::string, "String mask of zeros and ones, e.g. 000111000 indicating which coils to keep", "");
        NODE_PROPERTY(coils_out, int, "Number of coils to keep, coils with higher indices will be discarded", 128); // TODO: re-add limits
        std::vector<unsigned short> coil_mask_;
        unsigned int coils_in_;
        unsigned int coils_out_;
    };
}
