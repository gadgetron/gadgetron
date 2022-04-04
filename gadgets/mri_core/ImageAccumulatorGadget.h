/**
    \brief  
    \author Original: David Christoffer Hansen
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Untested 
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"
#include <mri_core_data.h>

namespace Gadgetron{
  class ImageAccumulatorGadget : public Core::ChannelGadget<IsmrmrdImageArray> 
    {
      public:
        using Core::ChannelGadget<IsmrmrdImageArray>::ChannelGadget;
        ImageAccumulatorGadget(const Core::Context& context, const Core::GadgetProperties& props);
        ~ImageAccumulatorGadget() override = default;
        void process(Core::InputChannel<IsmrmrdImageArray>& input, Core::OutputChannel& output) override;
      protected:
        NODE_PROPERTY(accumulate_dimension,std::string,"Dimension over which the images will be collected","contrast");
        NODE_PROPERTY(combine_along,std::string,"Dimension used for stacking the images","N");
      private:
        template<class T> auto extract_value(T& val);
        IsmrmrdImageArray combine_images(std::vector<IsmrmrdImageArray>&);
        size_t encoding_spaces;
        std::vector<IsmrmrdImageArray> images;
        std::vector<uint16_t> required_values;
        std::set<uint16_t> seen_values;
    };
}
