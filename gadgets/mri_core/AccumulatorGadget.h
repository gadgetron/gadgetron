/**
    \brief
    \test   Tested by: fs_csi.cfg
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "Node.h"

namespace Gadgetron{
  class AccumulatorGadget : public Core::ChannelGadget<mrd::Acquisition>
    {
      public:
        using Core::ChannelGadget<mrd::Acquisition>::ChannelGadget;
        AccumulatorGadget(const Core::Context& context, const Core::GadgetProperties& props);
        ~AccumulatorGadget() override;
        void process(Core::InputChannel<mrd::Acquisition>& input, Core::OutputChannel& output) override;
      protected:
        NODE_PROPERTY(image_series, int, "Image series", 0);
        hoNDArray<std::complex<float>>* buffer_;
        std::vector<size_t> dimensions_;
        float field_of_view_[3];
        size_t slices_;
        long long image_counter_;
        long long image_series_;
    };
}