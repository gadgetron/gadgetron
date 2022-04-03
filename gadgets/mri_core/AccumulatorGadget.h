/**
    \brief  
    \author Original: Thomas Sangild Sorensen
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Untested
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"

namespace Gadgetron{
  class AccumulatorGadget : public Core::ChannelGadget<Core::Acquisition> 
    {
      public:
        using Core::ChannelGadget<Core::Acquisition>::ChannelGadget;
        AccumulatorGadget(const Core::Context& context, const Core::GadgetProperties& props);
        ~AccumulatorGadget() override;
        void process(Core::InputChannel<Core::Acquisition>& input, Core::OutputChannel& output) override;
      protected:
        NODE_PROPERTY(image_series, int, "Image series", 0);
        hoNDArray<std::complex<float>>* buffer_;
        std::vector<size_t> dimensions_;
        std::vector<float> field_of_view_;
        size_t slices_;
        long long image_counter_;
        long long image_series_;
    };
}
