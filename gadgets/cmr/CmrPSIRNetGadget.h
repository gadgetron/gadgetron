/** \file   CmrPSIRNetGadget.h
    \brief  This gadget integrates the PSIRNet reconstruction for cardiac MRI.

            It is to process the PSIR LGE (Late Gadolinium Enhancement) imaging data. IR and PD are on dimension SET.

            S will be SET. N is average. 

            E2 is 1 for only supporing 2D imaging.

    \author Hui Xue
*/

#pragma once

#include "gadgetron_cmr_export.h"
#include "generic_recon_gadgets/GenericReconGadget.h"
#include "python_toolbox.h"

namespace Gadgetron {

    class EXPORTGADGETSCMR CmrPSIRNetGadget : public GenericReconGadget
    {
    public:
        GADGET_DECLARE(CmrPSIRNetGadget);

        typedef GenericReconGadget BaseClass;

        CmrPSIRNetGadget();
        ~CmrPSIRNetGadget();

        /// parameters for workflow
        GADGET_PROPERTY(model, std::string, "model file", "psirnet_model.pts");

        GADGET_PROPERTY(send_out_mag_IR, bool, "Whether to set out magIR images", true);
        GADGET_PROPERTY(offset_factor_after_SCC, double, "Offset factor after psir", 2048);
        GADGET_PROPERTY(scale_factor_after_SCC, double, "Scaling factor after psir", 1000);

    protected:

        bool prepare_AI();

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        std::vector<float> TI_;

        boost::python::object model_;
        bool model_loaded_;

        // gadgetron home
        std::string gt_home_;
        std::string model_dir_;

        // the raw recon results
        // [RO E1 E2 1 N S SLC]
        IsmrmrdImageArray res_psir_;
        IsmrmrdImageArray res_magir_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);

        // --------------------------------------------------
        // recon step functions
        // --------------------------------------------------
        virtual void perform_psir(IsmrmrdReconBit& recon_bit, size_t encoding);

        // compute window level for psir
        bool calculate_window_level(hoNDArray<std::complex<float>>& magPDFiltered, hoNDArray<std::complex<float>>& PSIRImage, float& window_center, float& window_width);

        // compute image header for PSIR and mag IR images
        int compute_image_header_psir_magir(IsmrmrdImageArray& res_psir, IsmrmrdImageArray& res_magir, size_t encoding);

        // --------------------------------------------------
        // overload functions
        // --------------------------------------------------
    };
}
