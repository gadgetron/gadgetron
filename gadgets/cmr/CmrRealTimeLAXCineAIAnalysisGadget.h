/** \file   CmrRealTimeLAXCineAIAnalysisGadget.h
    \brief  This is the class gadget for landmark detection of LAX RealTime Cine images.
    \author Hui Xue
*/

#pragma once

#include "gadgetron_cmr_export.h"
#include "GenericReconBase.h"
#include "hoMRImage.h"
#include "hoNDImageContainer2D.h"

namespace Gadgetron {

    class EXPORTGADGETSCMR CmrRealTimeLAXCineAIAnalysisGadget : public GenericReconImageBase
    {
    public:
        GADGET_DECLARE(CmrRealTimeLAXCineAIAnalysisGadget);

        typedef GenericReconImageBase BaseClass;
        typedef hoNDImageContainer2D < hoMRImage<float, 2> > ImageContainerMagType;

        CmrRealTimeLAXCineAIAnalysisGadget();
        ~CmrRealTimeLAXCineAIAnalysisGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the AI
        /// ------------------------------------------------------------------------------------
        GADGET_PROPERTY(perform_AI, bool, "Whether to perform AI", true);

        GADGET_PROPERTY(lax_landmark_detection_model, std::string, "Cine AI model file for lax landmark detection", "CMR_landmark_network_RO_352_E1_352_ch2_ch3_ch4_myo_pts_with_T1_LGE_LossMultiSoftProb_KLD_Dice_Pytorch_1.5.0_2020-06-17_20200617_111642.pts");
        GADGET_PROPERTY(oper_RO, size_t, "Operation image size for AI, RO", 352);
        GADGET_PROPERTY(oper_E1, size_t, "Operation image size for AI, E1", 352);
        GADGET_PROPERTY(pixel_size_send, double, "Pixel size used for AI and image sending", 1.0);

    protected:

        std::string gt_home_;
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
        virtual int perform_LAX_detection_AI(IsmrmrdImageArray& lax, IsmrmrdImageArray& lax_ai);

        // utility functions
        void convert_array_to_image_container(IsmrmrdImageArray& lax, hoNDImageContainer2D < hoMRImage<float, 2> >& lax_container);
        void convert_image_container_to_array(hoNDImageContainer2D < hoMRImage<float, 2> >& lax_container, IsmrmrdImageArray& lax);
        void plot_landmarks_on_images(hoNDImageContainer2D < hoMRImage<float, 2> > & lax_container, const hoNDArray<float>& pts);
        void attach_info_to_report(hoNDImageContainer2D < hoMRImage<float, 2> >& lax_container, const hoNDArray<float>& pts);
    };
}
