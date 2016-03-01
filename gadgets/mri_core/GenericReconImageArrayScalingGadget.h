/** \file   GenericReconGenericReconImageArrayScalingGadget.h
    \brief  This is the class gadget to scale the incoming IsmrmrdReconRes.

    \author Hui Xue
*/

#pragma once

#include <complex>
#include "gadgetron_mricore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"
#include "GadgetronTimer.h"

#include "GadgetStreamController.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_data.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconImageArrayScalingGadget : public Gadget1<IsmrmrdImageArray>
    {
    public:
        GADGET_DECLARE(GenericReconImageArrayScalingGadget);

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef Gadget1<IsmrmrdImageArray> BaseClass;

        GenericReconImageArrayScalingGadget();
        ~GenericReconImageArrayScalingGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------
        /// image scaling
        GADGET_PROPERTY(use_constant_scalingFactor, bool, "Whether to use constrant scaling; if not, the auto-scaling factor will be computed only ONCE", true);
        GADGET_PROPERTY(scalingFactor, float, "Default scaling ratio", 4.0);
        GADGET_PROPERTY(min_intensity_value, int, "Minimal intensity value for auto image scaling", 64);
        GADGET_PROPERTY(max_intensity_value, int, "Maxmimal intensity value for auto image scaling", 4095);
        GADGET_PROPERTY(auto_scaling_only_once, bool, "Whether to compute auto-scaling factor only once; if false, an auto-scaling factor is computed for every incoming image array", true);

        GADGET_PROPERTY(use_dedicated_scalingFactor_meta_field, std::string, "If this meta field exists, scale the images with the dedicated scaling factor", "Use_dedicated_scaling_factor");
        GADGET_PROPERTY(scalingFactor_dedicated, float, "Dedicated scaling ratio", 100.0);

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // number of encoding spaces in the protocol
        size_t num_encoding_spaces_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // number of times the process function is called
        size_t process_called_times_;

        // scaling factor used for every encoding space
        std::vector<double> scaling_factor_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_;

        // in verbose mode, more info is printed out
        bool verbose_;

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1);

        // scale the recon images
        virtual int compute_and_apply_scaling_factor(IsmrmrdImageArray& res, size_t encoding);

        // close call
        int close(unsigned long flags);
    };
}
