/** \file   GenericReconKSpaceFilteringGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian kspace filtering, working on the IsmrmrdImageArray.

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

#include "gtPlusIOAnalyze.h"

#include "GadgetStreamController.h"
#include "mri_core_data.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconKSpaceFilteringGadget : public Gadget1<IsmrmrdImageArray>
    {
    public:
        GADGET_DECLARE(GenericReconKSpaceFilteringGadget);

        typedef Gadget1<IsmrmrdImageArray> BaseClass;

        GenericReconKSpaceFilteringGadget();
        ~GenericReconKSpaceFilteringGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------
        GADGET_PROPERTY(skip_processing_meta_field, std::string, "If this meta field exists, pass the incoming image array to next gadget without processing", "Skip_processing_after_recon");

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

        /// ------------------------------------------------------------------------------------
        /// kspace filter parameters
        GADGET_PROPERTY_LIMITS(filterRO, std::string, "Kspace filter for RO dimension", "Gaussian",
            GadgetPropertyLimitsEnumeration, "Gaussian", "Hanning", "TaperedHanning", "None");

        GADGET_PROPERTY(filterRO_sigma, double, "Filter sigma for gaussian for RO dimension", 1.0);
        GADGET_PROPERTY(filterRO_width, double, "Filter width for tapered hanning for RO dimension", 0.15);

        // ------------------------------------------------------------------------------------

        GADGET_PROPERTY_LIMITS(filterE1, std::string, "Kspace filter for E1 dimension", "Gaussian",
            GadgetPropertyLimitsEnumeration, "Gaussian", "Hanning", "TaperedHanning", "None");

        GADGET_PROPERTY(filterE1_sigma, double, "Filter sigma for gaussian for E1 dimension", 1.0);
        GADGET_PROPERTY(filterE1_width, double, "Filter width for tapered hanning for E1 dimension", 0.15);

        // ------------------------------------------------------------------------------------

        GADGET_PROPERTY_LIMITS(filterE2, std::string, "Kspace filter for E2 dimension", "Gaussian",
            GadgetPropertyLimitsEnumeration, "Gaussian", "Hanning", "TaperedHanning", "None");

        GADGET_PROPERTY(filterE2_sigma, double, "Filter sigma for gaussian for E2 dimension", 1.0);
        GADGET_PROPERTY(filterE2_width, double, "Filter width for tapered hanning for E2 dimension", 0.15);

        // ------------------------------------------------------------------------------------

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // number of encoding spaces in the protocol
        size_t num_encoding_spaces_;

        // encoding FOV and recon FOV
        std::vector< std::vector<float> > encoding_FOV_;
        std::vector< std::vector<float> > recon_FOV_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // kspace filter for every encoding space
        std::vector< hoNDArray< std::complex<float> > > filter_RO_;
        std::vector< hoNDArray< std::complex<float> > > filter_E1_;
        std::vector< hoNDArray< std::complex<float> > > filter_E2_;

        // kspace buffer
        hoNDArray< std::complex<float> > kspace_buf_;

        // results of filtering
        hoNDArray< std::complex<float> > filter_res_;

        // number of times the process function is called
        size_t process_called_times_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // debug folder
        std::string debug_folder_full_path_;

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        // exporter
        Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

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

        // find kspace sampled range
        void find_kspace_sampled_range(size_t min, size_t max, size_t len, size_t& r);
    };
}
