/**
    \brief  Autoscales real-type images based on a given max value and the 99th percentile of the data
    \test   Tested by: epi_2d.cfg and others
*/

#pragma once

#include "PureGadget.h"
#include "Types.h"
#include "hoNDArray_math.h"
#include <algorithm>

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
