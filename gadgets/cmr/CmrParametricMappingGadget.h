/** \file   CmrParametricMappingGadget.h
    \brief  This is the class gadget for cardiac parametric mapping, working on the IsmrmrdImageArray.
    \author Hui Xue
*/

#pragma once

#include "gadgetron_cmr_export.h"
#include "generic_recon_gadgets/GenericReconBase.h"

namespace Gadgetron {

    class EXPORTGADGETSCMR CmrParametricMappingGadget : public GenericReconImageBase
    {
    public:
        GADGET_DECLARE(CmrParametricMappingGadget);

        typedef GenericReconImageBase BaseClass;

        CmrParametricMappingGadget();
        ~CmrParametricMappingGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the mapping
        /// ------------------------------------------------------------------------------------

        GADGET_PROPERTY(skip_processing_meta_field, std::string, "If this meta field exists, pass the incoming image array to next gadget without processing", GADGETRON_SKIP_PROCESSING_AFTER_RECON);

        GADGET_PROPERTY(imaging_prep_time_from_protocol, bool, "If true, read in imaging prep time from protocols; if false, read in from meta fields", true);

        // -------------------------------------------

        GADGET_PROPERTY(send_map, bool, "Whether to send out maps", true);
        GADGET_PROPERTY(send_sd_map, bool, "Whether to send out sd maps", true);

        GADGET_PROPERTY(color_lut_map, std::string, "Color lookup table for map", "GadgetronParametricMap.pal");
        GADGET_PROPERTY(window_center_map, double, "Window center for map", 4.0);
        GADGET_PROPERTY(window_width_map, double, "Window width for map", 8.0);

        GADGET_PROPERTY(color_lut_map_3T, std::string, "Color lookup table for map at 3T", "GadgetronParametricMap_3T.pal");
        GADGET_PROPERTY(window_center_map_3T, double, "Window center for map at 3T", 4.0);
        GADGET_PROPERTY(window_width_map_3T, double, "Window width for map at 3T", 8.0);

        GADGET_PROPERTY(scaling_factor_map, double, "Scale factor for map", 10.0);

        // -------------------------------------------

        GADGET_PROPERTY(color_lut_sd_map, std::string, "Color lookup table for sd map", "GadgetronParametricSDMap.pal");
        GADGET_PROPERTY(window_center_sd_map, double, "Window center for sd map", 4.0);
        GADGET_PROPERTY(window_width_sd_map, double, "Window width for sd map", 8.0);
        GADGET_PROPERTY(scaling_factor_sd_map, double, "Scale factor for sd map", 100.0);

        GADGET_PROPERTY(perform_hole_filling, bool, "Whether to perform hole filling on map", true);
        GADGET_PROPERTY(max_size_hole, int, "Maximal size for hole", 20);

        GADGET_PROPERTY(std_thres_masking, double, "Number of noise std for masking", 3.0);
        GADGET_PROPERTY(mapping_with_masking, bool, "Whether to compute and apply a mask for mapping", true);

        // ------------------------------------------------------------------------------------

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // field strength in T
        float field_strength_T_;

        // imaging prep time (e.g. inverison/saturation/echo time)
        std::vector<float> prep_times_;

        // encoding space size
        ISMRMRD::EncodingCounters meas_max_idx_;

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1);

        // close call
        int close(unsigned long flags);

        // function to perform the mapping
        // data: input image array [RO E1 E2 CHA N S SLC]
        // map and map_sd: mapping result and its sd
        // para and para_sd: other parameters of mapping and its sd
        virtual int perform_mapping(IsmrmrdImageArray& data, IsmrmrdImageArray& map, IsmrmrdImageArray& para, IsmrmrdImageArray& map_sd, IsmrmrdImageArray& para_sd) = 0;

        // fill image header and meta for maps
        virtual int fill_map_header(IsmrmrdImageArray& map);
        virtual int fill_sd_header(IsmrmrdImageArray& map_sd);

        // compute image mask
        virtual void compute_mask_for_mapping(const hoNDArray<float> &mag, hoNDArray<float> &mask,
                                              float scale_factor);
    };
}
