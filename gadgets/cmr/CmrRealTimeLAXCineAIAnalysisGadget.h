/** \file   CmrRealTimeLAXCineAIAnalysisGadget.h
    \brief  This is the class gadget for landmark detection of LAX RealTime Cine images.
    \author Hui Xue
*/

#pragma once

#include "generic_recon_gadgets/GenericReconBase.h"
#include "hoMRImage.h"
#include "hoNDImageContainer2D.h"
#include "python_toolbox.h"

namespace Gadgetron {

    class CmrRealTimeLAXCineAIAnalysisGadget : public GenericReconImageBase
    {
    public:
        typedef GenericReconImageBase BaseClass;
        typedef hoNDImageContainer2D < hoMRImage<float, 2> > ImageContainerMagType;

        CmrRealTimeLAXCineAIAnalysisGadget();
        ~CmrRealTimeLAXCineAIAnalysisGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the AI
        /// ------------------------------------------------------------------------------------
        GADGET_PROPERTY(perform_AI, bool, "Whether to perform AI", true);

        GADGET_PROPERTY(model_file, std::string, "Cine AI ONNX model file for lax landmark detection", "CMR_landmark_network_RO_320_E1_320_CH2_CH3_CH4_Myo_Pts_sFOV_LossMultiSoftProb_KLD_Dice_Pytorch_1.8.0a0+37c1f4a_2021-08-08_20210808_085042.onnx");
        GADGET_PROPERTY(model_file_sha256, std::string, "Checksum sha256 of the model to check its integrity", "48efe3e70b1ff083c9dd0066469f62bf495e52857d68893296e7375b69f691e4");
        GADGET_PROPERTY(oper_RO, size_t, "Operation image size for AI, RO", 320);
        GADGET_PROPERTY(oper_E1, size_t, "Operation image size for AI, E1", 320);
        GADGET_PROPERTY(pixel_size_send, double, "Pixel size used for AI and image sending", 1.0);

        GADGET_PROPERTY(model_url, std::string, "url to download the model", "https://gadgetrondata.blob.core.windows.net/cmr-ai-models/");
        GADGET_PROPERTY(model_dest, std::string, "destination folder to install model, under ${GADGETRON_INSTALL_PYTHON_MODULE_PATH}", "cmr_lax_landmark_detection");

    protected:

        std::string gt_home_;
        std::string gt_model_home_;

        mrd::EncodingCounters meas_max_idx_;

        boost::python::object model_;

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(const mrd::Header& header);
        virtual int process(Gadgetron::GadgetContainerMessage< mrd::ImageArray >* m1);

        // close call
        int close(unsigned long flags);

        // function to perform the mapping
        virtual int perform_LAX_detection_AI(mrd::ImageArray& lax, mrd::ImageArray& lax_ai);

        // utility functions
        void convert_array_to_image_container(mrd::ImageArray& lax, hoNDImageContainer2D < hoMRImage<float, 2> >& lax_container);
        void convert_image_container_to_array(hoNDImageContainer2D < hoMRImage<float, 2> >& lax_container, mrd::ImageArray& lax);
        void plot_landmarks_on_images(hoNDImageContainer2D < hoMRImage<float, 2> > & lax_container, const hoNDArray<float>& pts);
        void attach_info_to_report(hoNDImageContainer2D < hoMRImage<float, 2> >& lax_container, const hoNDArray<float>& pts);
    };
}
