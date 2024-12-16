#pragma once

#include "PureGadget.h"
namespace Gadgetron {
    using PercentileScaleImageTypes = std::variant<mrd::Image<float>, mrd::ImageArray>;
/***
 * This Gadget rescales magnitude images so that their 99% percentile becomes max_value
 */
    class ScaleGadget : public Core::PureGadget<PercentileScaleImageTypes, PercentileScaleImageTypes > {
    public:
        using PureGadget<PercentileScaleImageTypes, PercentileScaleImageTypes>::PureGadget;
       PercentileScaleImageTypes process_function(PercentileScaleImageTypes args) const override;

    protected:
        NODE_PROPERTY(max_value, float, "Percentile value (after scaling)", 2048);
        NODE_PROPERTY(percentile,float,"Percentile to use.",98.5);
    };
}


