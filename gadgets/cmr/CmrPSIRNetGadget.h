/** \file   CmrPSIRNetGadget.h
    \brief  This gadget integrates the PSIRNet reconstruction for cardiac MRI.

            It is to process the PSIR LGE (Late Gadolinium Enhancement) imaging data. IR and PD are on dimension SET.

            S will be SET. [RO, E1, E2, CHA, N, S, SLC]

            E2 is 1 for only supporing 2D imaging.

    \author Hui Xue
*/

#pragma once

#include "gadgetron_cmr_export.h"
#include "generic_recon_gadgets/GenericReconGadget.h"

namespace Gadgetron {

    class EXPORTGADGETSCMR CmrPSIRNetGadget : public GenericReconGadget
    {
    public:
        GADGET_DECLARE(CmrPSIRNetGadget);

        typedef GenericReconGadget BaseClass;

        CmrPSIRNetGadget();
        ~CmrPSIRNetGadget();

        /// parameters for workflow
        GADGET_PROPERTY(send_out_mag_IR, bool, "Whether to set out magIR images", true);
        GADGET_PROPERTY(model, std::string, "model file", "psirnet_model.pts");

    protected:

        bool prepare_AI();

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        boost::python::object model_;
        bool model_loaded_;

        // gadgetron home
        std::string gt_home_;
        std::string model_dir_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);

        // --------------------------------------------------
        // recon step functions
        // --------------------------------------------------
        virtual void perform_recon(IsmrmrdReconBit& recon_bit, size_t encoding);

        // --------------------------------------------------
        // overload functions
        // --------------------------------------------------
        // send out the recon results
        virtual int prep_image_header_send_out(IsmrmrdImageArray& res, size_t n, size_t s, size_t slc, size_t encoding, int series_num, const std::string& data_role);
    };
}
