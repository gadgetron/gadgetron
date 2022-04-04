/**
    \brief  Rescales magnitude images so that their 99% percentile becomes max_value
    \author Original: David Christopher Hansen
    \test   Untested
*/

#pragma once

#include "PureGadget.h"
#include <mri_core_data.h>
namespace Gadgetron {

    using PercentileScaleImageTypes = Core::variant<Core::Image<float>, IsmrmrdImageArray>;
    class ScaleGadget : public Core::PureGadget<PercentileScaleImageTypes, PercentileScaleImageTypes > {
        
    public:
        using PureGadget<PercentileScaleImageTypes, PercentileScaleImageTypes>::PureGadget;
       PercentileScaleImageTypes process_function(PercentileScaleImageTypes args) const override;

    protected:
        NODE_PROPERTY(max_value, float, "Percentile value (after scaling)", 2048);
        NODE_PROPERTY(percentile,float,"Percentile to use.",98.5);

    };
}


