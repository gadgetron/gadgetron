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

     mrd::ImageArray autoscale(mrd::ImageArray images, float max_value, float percentile) {
         percentile_scale(images.data, max_value, percentile);
         return images;
    }

    mrd::Image<float> autoscale(mrd::Image<float> image, float max_value, float percentile) {
        if (image.head.image_type != mrd::ImageType::kMagnitude) {
            return image;
        }

        percentile_scale(image.data, max_value, percentile);
        return image;
    }

} // namespace

GADGETRON_GADGET_EXPORT(ScaleGadget)

PercentileScaleImageTypes ScaleGadget::process_function(PercentileScaleImageTypes args) const {
    return std::visit([&](auto image) { return PercentileScaleImageTypes(autoscale(image, max_value, percentile)); },
                      args);
}

} // namespace Gadgetron
