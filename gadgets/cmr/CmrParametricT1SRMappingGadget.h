/**
\file   CmrParametricT1SRMappingGadget.h
\brief  This is the class gadget for cardiac T1 SR mapping, working on the IsmrmrdImageArray.
\author Hui Xue
*/

#pragma once

#include "CmrParametricMappingGadget.h"

namespace Gadgetron {

    class EXPORTGADGETSCMR CmrParametricT1SRMappingGadget : public CmrParametricMappingGadget
    {
    public:
        GADGET_DECLARE(CmrParametricT1SRMappingGadget);

        typedef CmrParametricMappingGadget BaseClass;

        CmrParametricT1SRMappingGadget();
        ~CmrParametricT1SRMappingGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the mapping
        /// ------------------------------------------------------------------------------------

        GADGET_PROPERTY(max_iter, size_t, "Maximal number of iterations", 150);
        GADGET_PROPERTY(thres_func, double, "Threshold for minimal change of cost function", 1e-4);
        GADGET_PROPERTY(max_T1, double, "Maximal T1 allowed in mapping (ms)", 4000);

        GADGET_PROPERTY(anchor_image_index, size_t, "Index for anchor image; by default, the first image is the anchor (without SR pulse)", 0);
        GADGET_PROPERTY(anchor_TS, double, "Saturation time for anchor", 10000);

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);

        virtual int close(unsigned long flags);

        // function to perform the mapping
        // data: input image array [RO E1 E2 CHA N S SLC]
        // map and map_sd: mapping result and its sd
        // para and para_sd: other parameters of mapping and its sd
        virtual int perform_mapping(IsmrmrdImageArray& data, IsmrmrdImageArray& map, IsmrmrdImageArray& para, IsmrmrdImageArray& map_sd, IsmrmrdImageArray& para_sd);
    };
}
