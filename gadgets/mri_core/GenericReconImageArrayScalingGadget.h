/** \file   GenericReconGenericReconImageArrayScalingGadget.h
    \brief  This is the class gadget to scale the incoming IsmrmrdReconRes.

    \author Hui Xue
*/

#pragma once

#include "GenericReconBase.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconImageArrayScalingGadget : public GenericReconImageBase
    {
    public:
        GADGET_DECLARE(GenericReconImageArrayScalingGadget);

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef GenericReconImageBase BaseClass;

        GenericReconImageArrayScalingGadget();
        ~GenericReconImageArrayScalingGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------
        /// image scaling
        GADGET_PROPERTY(use_constant_scalingFactor, bool, "Whether to use constraint scaling; if not, the auto-scaling factor will be computed only ONCE", true);
        GADGET_PROPERTY(scalingFactor, float, "Default scaling ratio", 4.0);
        GADGET_PROPERTY(min_intensity_value, int, "Minimal intensity value for auto image scaling", 64);
        GADGET_PROPERTY(max_intensity_value, int, "Maximal intensity value for auto image scaling", 4095);
        GADGET_PROPERTY(auto_scaling_only_once, bool, "Whether to compute auto-scaling factor only once; if false, an auto-scaling factor is computed for every incoming image array", true);

        GADGET_PROPERTY(use_dedicated_scalingFactor_meta_field, std::string, "If this meta field exists, scale the images with the dedicated scaling factor", "Use_dedicated_scaling_factor");
        GADGET_PROPERTY(scalingFactor_dedicated, float, "Dedicated scaling ratio", 100.0);
        GADGET_PROPERTY(scalingFactor_gfactor_map, float, "Scaling ratio for gfactor map", 100.0);
        GADGET_PROPERTY(scalingFactor_snr_map, float, "Scaling ratio for snr map", 10.0);
        GADGET_PROPERTY(scalingFactor_snr_std_map, float, "Scaling ratio for snr standard deviation map", 1000.0);

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // scaling factor used for every encoding space
        std::vector<double> scaling_factor_;

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1);

        // scale the recon images
        virtual int compute_and_apply_scaling_factor(IsmrmrdImageArray& res, size_t encoding);

    };
}
