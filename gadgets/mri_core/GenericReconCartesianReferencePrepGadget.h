/** \file   GenericReconCartesianReferencePrepGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian reconstruction to prepare the reference data, working on the IsmrmrdReconData.
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

#include "hoNDArray_utils.h"

//#include "gtPlusIOAnalyze.h"

#include "GadgetStreamController.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"

#include "mri_core_data.h"
#include "mri_core_utility.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconCartesianReferencePrepGadget : public Gadget1<IsmrmrdReconData>
    {
    public:
        GADGET_DECLARE(GenericReconCartesianReferencePrepGadget);

        typedef Gadget1<IsmrmrdReconData> BaseClass;

        GenericReconCartesianReferencePrepGadget();
        ~GenericReconCartesianReferencePrepGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

        /// ref preparation
        /// whether to average all N for ref generation
        /// for the interleaved mode, the sampling times will be counted and used for averaging
        /// it is recommended to set N as the interleaved dimension
        GADGET_PROPERTY(average_all_ref_N, bool, "Whether to average all N for ref generation", true);
        /// whether to average all S for ref generation
        GADGET_PROPERTY(average_all_ref_S, bool, "Whether to average all S for ref generation", false);
        /// whether to update ref for every incoming IsmrmrdReconData; for some applications, we may want to only compute ref data once
        /// if false, the ref will only be prepared for the first incoming IsmrmrdReconData
        GADGET_PROPERTY(prepare_ref_always, bool, "Whether to prepare ref for every incoming IsmrmrdReconData", true);

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // number of encoding spaces in the protocol
        size_t num_encoding_spaces_;
        /// indicate whether ref has been prepared for an encoding space
        std::vector<bool> ref_prepared_;

        // for every encoding space
        // calibration mode
        std::vector<Gadgetron::ismrmrdCALIBMODE> calib_mode_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // number of times the process function is called
        size_t process_called_times_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // in verbose mode, more info is printed out
        bool verbose_;

        //// debug folder
        //std::string debug_folder_full_path_;

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        //// exporter
        //Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);
    };
}
