/** \file   FatWaterSeparationGadget.h
    \brief  The gadget performs the fat water seperation, calling the toolbox functionalities
    \author Hui Xue
*/

#pragma once

#include "GenericReconImagePostProcessingGadget.h"

namespace Gadgetron { 

class EXPORTGADGETSMRICORE FatWaterSeparationGadget : public GenericReconImagePostProcessing3DGadget
{
public:
    GADGET_DECLARE(FatWaterSeparationGadget);

    typedef GenericReconImagePostProcessing3DGadget BaseClass;

    typedef BaseClass::ValueType ValueType;
    typedef BaseClass::ImageType ImageType;
    typedef BaseClass::ImageMagType ImageMagType;
    typedef BaseClass::ImageBufferType ImageBufferType;
    typedef BaseClass::ImgArrayType ImgArrayType;

    FatWaterSeparationGadget();
    ~FatWaterSeparationGadget();

    virtual int close(unsigned long flags);

protected:

    GADGET_PROPERTY(send_multi_echo_images, bool, "Whether to send out multi-echo images", false);
    GADGET_PROPERTY(send_water_images, bool, "Whether to send out water images", false);
    GADGET_PROPERTY(send_fat_images, bool, "Whether to send out fat images", false);
    GADGET_PROPERTY(send_t2_star_map, bool, "Whether to send out t2* map", false);
    GADGET_PROPERTY(send_field_map, bool, "Whether to send out field map", false);

    GADGET_PROPERTY(t2_star_color_map, std::string, "Color map for t2* map", "GadgetronColorMap.pal");
    GADGET_PROPERTY(t2_star_scaling_factor, double, "Scaling factor for t2* map", 100.0);
    GADGET_PROPERTY(t2_star_window_center, double, "Window center in ms", 25.0);
    GADGET_PROPERTY(t2_star_window_width, double, "Window width in ms", 50.0);

    GADGET_PROPERTY(field_map_color_map, std::string, "Color map for field map", "GadgetronColorMap.pal");
    GADGET_PROPERTY(field_map_scaling_factor, double, "Scaling factor for field map map", 10.0);
    GADGET_PROPERTY(field_map_offset, double, "Offset for field map map", 4096);
    GADGET_PROPERTY(field_map_window_center, double, "Window center in Hz", 25.0);
    GADGET_PROPERTY(field_map_window_width, double, "Window width in Hz", 50.0);

    virtual int process_config(ACE_Message_Block* mb);
    virtual int process_image_buffer(ImageBufferType& ori);

    /// perform the fat water separation
    int perform_fat_water(ImageBufferType& input, ImageBufferType& water, ImageBufferType& fat, ImageBufferType& t2s_map, ImageBufferType& field_map);
    int perform_fat_water(ImageContainerType& input, ImageContainerType& water, ImageContainerType& fat, ImageContainerType& t2s_map, ImageContainerType& field_map);

    // echo times in ms
    std::vector<double> echo_times_;
};

}
