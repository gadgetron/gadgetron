/**
    \brief  Autoscales real-type images based on a given max value and the 99th percentile of the data
    \author Original: David Christoffer Hansen
    \author PureGadget Conversion: Andrew Dupuis
    \test   Untested
*/

#pragma once

#include "PureGadget.h"
#include "Types.h"

namespace Gadgetron{
    class AutoScaleGadget : public Core::PureGadget<Core::AnyImage, Core::AnyImage> {
    public:
      using Core::PureGadget<Core::AnyImage,Core::AnyImage>::PureGadget;
        Core::AnyImage process_function(Core::AnyImage image) const override;
    protected:
        NODE_PROPERTY(max_value, float, "Maximum value (after scaling)", 2048);
        NODE_PROPERTY(histogram_bins, unsigned int, "Number of Histogram Bins", 100);
        float current_scale_;
        std::vector<size_t> histogram_;
    };
}
