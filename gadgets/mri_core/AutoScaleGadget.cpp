#include "AutoScaleGadget.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"

namespace Gadgetron {


    Core::Image<float> AutoScaleGadget::process_function(Core::Image<float> image) const {
        auto& header = std::get<ISMRMRD::ImageHeader>(image);
        auto& data   = std::get<hoNDArray<float>>(image);
        if (header.image_type != ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE)
            return image;

        auto scale = max_value / Gadgetron::percentile_approx(data,percentile/100);
        data *= scale;

        return image;
    }
    GADGETRON_GADGET_EXPORT(AutoScaleGadget)

}
