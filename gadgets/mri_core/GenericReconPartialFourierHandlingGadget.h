/** \file   GenericReconPartialFourierHandlingGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian partial fourier handling, working on the IsmrmrdImageArray.

            This class is a part of general cartesian recon chain.

            The sampled kspace region is found from image meta in fields:

            sampling_limits_RO
            sampling_limits_E1
            sampling_limits_E2

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

//#include "gtPlusIOAnalyze.h"

#include "GadgetStreamController.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDFFT.h"

#include "mri_core_data.h"
#include "mri_core_utility.h"
#include "mri_core_partial_fourier.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconPartialFourierHandlingGadget : public Gadget1<IsmrmrdImageArray>
    {
    public:
        GADGET_DECLARE(GenericReconPartialFourierHandlingGadget);

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef ValueType T;

        typedef Gadget1<IsmrmrdImageArray> BaseClass;

        GenericReconPartialFourierHandlingGadget();
        ~GenericReconPartialFourierHandlingGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------
        GADGET_PROPERTY(skip_processing_meta_field, std::string, "If this meta field exists, pass the incoming image array to next gadget without processing", "Skip_processing_after_recon");

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // number of encoding spaces in the protocol
        size_t num_encoding_spaces_;

        // acceleration factor for E1 and E2
        std::vector<double> acceFactorE1_;
        std::vector<double> acceFactorE2_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // sampled range
        size_t startRO_;
        size_t endRO_;
        size_t startE1_;
        size_t endE1_;
        size_t startE2_;
        size_t endE2_;

        // kspace buffer
        hoNDArray<T> kspace_buf_;

        // results of pf handling
        hoNDArray<T> pf_res_;

        // number of times the process function is called
        size_t process_called_times_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // debug folder
        //std::string debug_folder_full_path_;

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        // exporter
        //Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

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

        // --------------------------------------------------
        // implementation functions
        // --------------------------------------------------
        virtual int perform_partial_fourier_handling() = 0;

    };
}
