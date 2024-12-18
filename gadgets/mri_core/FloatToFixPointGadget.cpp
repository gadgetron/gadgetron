/*
*       FloatToFixPointGadget.cpp
*
*       Created on: March 10, 2014
*       Author: Hui Xue
*/

#include "FloatToFixPointGadget.h"
#include "mri_core_def.h"
#include "hoNDArray_math.h"


#include <boost/math/constants/constants.hpp>

namespace Gadgetron
{

    template<typename T, typename Base >
    void FloatToFixPointGadget<T, Base>::process(Core::InputChannel<mrd::Image<float>> &input, Core::OutputChannel &output) {

        auto self = static_cast<Base&>(*this);

        auto clamp = [&](float val){
            return lround(std::min<float>(std::max<float>(val,self.min_intensity),self.max_intensity));
        };
        auto magnitude = [&](auto val){
            return T(clamp(std::abs(val)));
        };

        auto real_value = [&](auto val){
            return T(clamp(float(val)+self.intensity_offset));
        };

        auto phase = [&](float val){
            return T(clamp((val*self.intensity_offset/boost::math::float_constants::pi)+self.intensity_offset));
        };

        for (auto image_in: input) {
            mrd::Image<T> image_out;
            image_out.head = image_in.head;
            image_out.meta = image_in.meta;
            image_out.data.create(image_in.data.dimensions());
            auto& output_data = image_out.data;

            // Now we're ready to transform the image data
            switch (image_in.head.image_type) {
                case mrd::ImageType::kMagnitude: {
                    std::transform(image_in.data.begin(), image_in.data.end(), output_data.begin(), magnitude);
                }
                    break;

                case mrd::ImageType::kReal:
                case mrd::ImageType::kImag: {
                    std::transform(image_in.data.begin(), image_in.data.end(), output_data.begin(), real_value);

                    if (image_out.meta.count(GADGETRON_IMAGE_WINDOWCENTER) && image_out.meta[GADGETRON_IMAGE_WINDOWCENTER].size() > 0) {
                        long windowCenter = std::get<long>(image_out.meta[GADGETRON_IMAGE_WINDOWCENTER].front());
                        image_out.meta[GADGETRON_IMAGE_WINDOWCENTER] = {windowCenter + self.intensity_offset};
                    }
                }
                    break;

                case mrd::ImageType::kPhase: {
                    std::transform(image_in.data.begin(), image_in.data.end(), output_data.begin(), phase);
                }
                    break;

                default:
                    throw std::runtime_error("Unknown image type in Image");
            }

            output.push(std::move(image_out));
        }
    }

    GADGETRON_GADGET_EXPORT(FloatToUShortGadget)
    GADGETRON_GADGET_EXPORT(FloatToShortGadget)
    GADGETRON_GADGET_EXPORT(FloatToIntGadget)
    GADGETRON_GADGET_EXPORT(FloatToUIntGadget)
}
