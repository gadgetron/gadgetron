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

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

        /// ------------------------------------------------------------------------------------
        /// partial fourier and asymmetric echo handling
        GADGET_PROPERTY_LIMITS(partial_fourier_algo, std::string, "Partial fourier handling method", "ZeroFillingFilter",
            GadgetPropertyLimitsEnumeration, "None", "ZeroFillingFilter", "POCS");

        // ------------------------------------------------------------------------------------

        GADGET_PROPERTY(partial_fourier_filter_RO_width, double, "Partial fourier filter width for tapered hanning for RO dimension", 0.15);
        GADGET_PROPERTY(partial_fourier_filter_E1_width, double, "Partial fourier filter width for tapered hanning for E1 dimension", 0.15);
        GADGET_PROPERTY(partial_fourier_filter_E2_width, double, "Partial fourier filter width for tapered hanning for E2 dimension", 0.15);

        GADGET_PROPERTY(partial_fourier_filter_densityComp, bool, "Whether to apply density compensation for RO dimension", false);

        // ------------------------------------------------------------------------------------

        GADGET_PROPERTY(partial_fourier_POCS_iters, int, "Number of iterations for POCS PF handling", 6);
        GADGET_PROPERTY(partial_fourier_POCS_thres, double, "Threshold for POSC PF handling", 0.01);
        GADGET_PROPERTY(partial_fourier_POCS_transitBand, int, "Transition band width for POCS PF handling", 24);
        GADGET_PROPERTY(partial_fourier_POCS_transitBand_E2, int, "Transition band width for POCS PF handling for E2 dimension", 16);

        // ------------------------------------------------------------------------------------

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

        // kspace buffer
        hoNDArray<T> kspace_buf_;

        // results of pf handling
        hoNDArray<T> pf_res_;

        // number of times the process function is called
        size_t process_called_times_;

        // partial fourier filters, avoid recomputing
        hoNDArray<T> filter_pf_RO_;
        hoNDArray<T> filter_pf_E1_;
        hoNDArray<T> filter_pf_E2_;

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
    };
}
