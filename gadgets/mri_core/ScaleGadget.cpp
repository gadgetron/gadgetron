#include "ScaleGadget.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include <mri_core_def.h>

namespace Gadgetron {

namespace {

    void percentile_scale(hoNDArray<std::complex<float>>& data, float max_value, float percentile){
        auto abs_image = abs(data);
        auto scale = max_value / Gadgetron::percentile(abs_image, percentile / 100);
        data *= scale;
    }
    void percentile_scale(hoNDArray<float>& data, float max_value, float percentile){
        auto scale = max_value / Gadgetron::percentile(data, percentile / 100);
        data *= scale;
    }

     IsmrmrdImageArray autoscale(IsmrmrdImageArray images, float max_value, float percentile) {
        percentile_scale(images.data_,max_value,percentile);
        return images;
    }
        Core::Image<float> autoscale(Core::Image<float> image,float max_value, float percentile) {
            auto& header = std::get<ISMRMRD::ImageHeader>(image);
            auto& data   = std::get<hoNDArray<float>>(image);
            if (header.image_type != ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE)
                return image;

            percentile_scale(data,max_value,percentile);
            return image;
        }
    }
    GADGETRON_GADGET_EXPORT(ScaleGadget)

    PercentileScaleImageTypes ScaleGadget::process_function(PercentileScaleImageTypes args) const {
        return Core::visit([&](auto image){ return PercentileScaleImageTypes(autoscale(image,max_value,percentile));},args);
    }
}
