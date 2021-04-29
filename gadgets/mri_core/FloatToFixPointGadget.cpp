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
    template<class T>
    constexpr auto ismrmrd_image_type(){
        if constexpr (std::is_same_v<T,unsigned short>) return ISMRMRD::ISMRMRD_USHORT;
        if constexpr (std::is_same_v<T,short>) return ISMRMRD::ISMRMRD_SHORT;
        if constexpr (std::is_same_v<T,int >) return ISMRMRD::ISMRMRD_INT;
        if constexpr (std::is_same_v<T,unsigned int>) return ISMRMRD::ISMRMRD_UINT;

        throw std::runtime_error("Unsupported type");

        }

    template<typename T, typename Base >
    void FloatToFixPointGadget<T,Base >::process(Core::InputChannel<Core::Image<float>> &input, Core::OutputChannel &output) {

        auto self = static_cast<Base&>(*this);

        auto clamp = [&](float val){
            return std::min<float>(std::max<float>(val,self.min_intensity),self.max_intensity);
        };
        auto magnitude = [&](auto val){
            return clamp(std::abs(val));
        };

        auto phase = [&](float val){
            return clamp((val*self.intensity_offset/boost::math::float_constants::pi)+self.intensity_offset);
        };


        for (auto [img_header,data,meta] : input) {
            GDEBUG("Float to shorts norm: %f \n", nrm2(&data));

            switch (img_header.image_type) {
                case ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE: {
                    std::transform(data.begin(),data.end(),data.begin(),magnitude);
                }
                    break;

                case ISMRMRD::ISMRMRD_IMTYPE_REAL:
                case ISMRMRD::ISMRMRD_IMTYPE_IMAG: {
                    std::transform(data.begin(),data.end(),data.begin(),clamp);

                    if (meta) {
                        if (meta->length(GADGETRON_IMAGE_WINDOWCENTER) > 0) {
                            long windowCenter;
                            windowCenter = meta->as_long(GADGETRON_IMAGE_WINDOWCENTER, 0);
                            meta->set(GADGETRON_IMAGE_WINDOWCENTER,
                                                    windowCenter + (long) self.intensity_offset);
                        }
                    }
                }
                    break;

                case ISMRMRD::ISMRMRD_IMTYPE_PHASE: {
                    std::transform(data.begin(),data.end(),data.begin(),phase );

                }
                    break;

                default:
                    throw std::runtime_error("Unknown image type in Image");

            }

            img_header.data_type = ismrmrd_image_type<T>();
            output.push(img_header,std::move(data),std::move(meta));
        }

    }


    GADGETRON_GADGET_EXPORT(FloatToUShortGadget)
    GADGETRON_GADGET_EXPORT(FloatToShortGadget)
    GADGETRON_GADGET_EXPORT(FloatToIntGadget)
    GADGETRON_GADGET_EXPORT(FloatToUIntGadget)
}
