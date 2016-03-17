/** \file   GenericReconNoiseStdMapComputingGadget.h
    \brief  This is the class gadget to compute standard deviation map, working on the IsmrmrdImageArray.

            This class is a part of general cartesian recon chain. It computes the std map on incoming SNR images.

\author     Hui Xue
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

// #include "gtPlusIOAnalyze.h"

#include "GadgetStreamController.h"

#include "mri_core_data.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconNoiseStdMapComputingGadget : public Gadget1<IsmrmrdImageArray>
    {
    public:
        GADGET_DECLARE(GenericReconNoiseStdMapComputingGadget);

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef Gadget1<IsmrmrdImageArray> BaseClass;

        GenericReconNoiseStdMapComputingGadget();
        ~GenericReconNoiseStdMapComputingGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

        /// ------------------------------------------------------------------------------------
        /// start N index to compute std map
        GADGET_PROPERTY(start_N_for_std_map, int, "Start N index to compute std map", 5);

        // ------------------------------------------------------------------------------------

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // number of encoding spaces in the protocol
        size_t num_encoding_spaces_;

        // number of times the process function is called
        size_t process_called_times_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // debug folder
        // std::string debug_folder_full_path_;

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        // exporter
        // Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

        // in verbose mode, more info is printed out
        bool verbose_;

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1);

        // close call
        int close(unsigned long flags);
    };
}
