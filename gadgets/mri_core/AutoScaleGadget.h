#pragma once

#include "PureGadget.h"
namespace Gadgetron {

    class AutoScaleGadget : public Core::TypedPureGadget<Core::Image<float>, Core::Image<float>> {
    public:
        using TypedPureGadget<Core::Image<float>,Core::Image<float>>::TypedPureGadget;
        Core::Image<float> process_function(Core::Image<float> args) const override;

    protected:
        NODE_PROPERTY(max_value, float, "Maximum value (after scaling)", 2048);
    };
}


