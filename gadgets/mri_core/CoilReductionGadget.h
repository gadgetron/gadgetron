/**
    \brief  Reduces number of coils in an acquisition based on a mask or threshold count
    \test   Tested by: simple_gre_3d.cfg, gpu_spiral_realtime_deblurring.cfg, gpu_fixed_radial_mode1_realtime.cfg, and others
*/

#pragma once

#include "Node.h"
#include "hoNDArray.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

namespace Gadgetron{
  class CoilReductionGadget : public Core::ChannelGadget<mrd::Acquisition>
    {
      public:
        using Core::ChannelGadget<mrd::Acquisition>::ChannelGadget;
        CoilReductionGadget(const Core::Context& context, const Core::GadgetProperties& props);
        ~CoilReductionGadget() override = default;
        void process(Core::InputChannel<mrd::Acquisition>& input, Core::OutputChannel& output) override;
      protected:
        NODE_PROPERTY(coil_mask, std::string, "String mask of zeros and ones, e.g. 000111000 indicating which coils to keep", "");
        NODE_PROPERTY(coils_out, int, "Number of coils to keep, coils with higher indices will be discarded", 128);
        std::vector<unsigned short> coil_mask_;
        unsigned int coils_in_;
        unsigned int coils_out_;
    };
}
