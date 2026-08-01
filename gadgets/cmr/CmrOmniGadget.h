/** \file   CmrOmniGadget.h
    \brief  This gadget integrates the OmniCMR reconstruction for cardiac MRI.

            if N > 1, recon will be long N dimension.
            If S > 1, N==1, recon will  be along S dimension.
            If S == 1 and N == 1, recon will be perform on each 2D image.

            E2 is 1 for only supporing 2D imaging.

    \author Hui Xue
*/

#pragma once

#include "gadgetron_cmr_export.h"
#include "generic_recon_gadgets/GenericReconGadget.h"
#include "python_toolbox.h"

namespace Gadgetron {

    class EXPORTGADGETSCMR CmrOmniGadget : public GenericReconGadget
    {
    public:
        GADGET_DECLARE(CmrOmniGadget);

        typedef GenericReconGadget BaseClass;

        CmrOmniGadget();
        ~CmrOmniGadget();

        /// parameters for workflow
        GADGET_PROPERTY(model, std::string, "model file", "omnicmr.pts");
        GADGET_PROPERTY(file_type, std::string, "file type, must be Perfusion, Retro, RT", "Perfusion");

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

        // the raw recon results
        // [RO E1 E2 1 N S SLC]
        IsmrmrdImageArray res_;

        size_t num_of_PD_images_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);

        // --------------------------------------------------
        // recon step functions
        // --------------------------------------------------
        virtual void perform_omni_recon(IsmrmrdReconBit& recon_bit, size_t encoding);

        // --------------------------------------------------
        // overload functions
        // --------------------------------------------------
    };
}