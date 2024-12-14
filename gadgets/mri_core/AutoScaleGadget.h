/**
    \brief  Autoscales real-type images based on a given max value and the 99th percentile of the data
    \test   Tested by: epi_2d.cfg and others
*/

#pragma once

#include "PureGadget.h"
#include "hoNDArray_math.h"
#include <algorithm>

namespace Gadgetron{
    class AutoScaleGadget : public Core::PureGadget<mrd::AnyImage, mrd::AnyImage> {
    public:
        using Core::PureGadget<mrd::AnyImage, mrd::AnyImage>::PureGadget;
        mrd::AnyImage process_function(mrd::AnyImage image) const override;
    protected:
        NODE_PROPERTY(max_value, float, "Maximum value (after scaling)", 2048);
        NODE_PROPERTY(histogram_bins, unsigned int, "Number of Histogram Bins", 100);
        float current_scale_;
        std::vector<size_t> histogram_;
    };
}
