#pragma once

#include "PureGadget.h"
namespace Gadgetron {
/***
 * This Gadget rescales magnitude images so that their 99% percentile becomes max_value
 */
    class AutoScaleGadget : public Core::PureGadget<Core::Image<float>, Core::Image<float>> {
    public:
        using PureGadget<Core::Image<float>,Core::Image<float>>::PureGadget;
        Core::Image<float> process_function(Core::Image<float> args) const override;

    protected:
        NODE_PROPERTY(max_value, float, "Percentile value (after scaling)", 2048);
        NODE_PROPERTY(percentile,float,"Percentile to use.",98.5);
    };
}


