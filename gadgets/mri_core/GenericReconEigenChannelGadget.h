/** \file   GenericReconEigenChannelGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian reconstruction to convert the data into eigen channel, working on the IsmrmrdReconData.
            If incoming data has the ref, ref data will be used to compute KLT coefficients
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

// #include "gtPlusIOAnalyze.h"

#include "GadgetStreamController.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDKLT.h"

#include "mri_core_data.h"
#include "mri_core_utility.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconEigenChannelGadget : public Gadget1<IsmrmrdReconData>
    {
    public:
        GADGET_DECLARE(GenericReconEigenChannelGadget);

        typedef Gadget1<IsmrmrdReconData> BaseClass;
        typedef hoNDKLT< std::complex<float> > KLTType;

        GenericReconEigenChannelGadget();
        ~GenericReconEigenChannelGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// compute KLT coefficients
        /// whether to average all N for coefficient computation
        /// for the interleaved mode, the sampling times will be counted and used for averaging
        GADGET_PROPERTY(average_all_ref_N, bool, "Whether to average all N for ref generation", true);
        /// whether to average all S for coefficient computation
        GADGET_PROPERTY(average_all_ref_S, bool, "Whether to average all S for ref generation", false);

        /// if update_eigen_channel_coefficients==true, every incoming IsmrmrdReconData will be used to compute KLT coefficients
        /// and the older one will be replaced
        /// if update_eigen_channel_coefficients==false, the KLT coefficients will be computed only once for the first incoming IsmrmrdReconData
        GADGET_PROPERTY(update_eigen_channel_coefficients, bool, "Whether to update KLT coeffients for eigen channel computation", false);

        /// optionally, upstream coil compression can be applied
        /// if upstream_coil_compression==true, only kept channels will be sent out to next gadgets and other channels will be removed
        /// no matter whether upstream_coil_compression is true or false, all channels will be converted into eigen channel 
        GADGET_PROPERTY(upstream_coil_compression, bool, "Whether to perform upstream coil compression", true);
        /// the logic here is that if upstream_coil_compression_num_modesKept>0, only upstream_coil_compression_num_modesKept channels will be kept
        /// if upstream_coil_compression_num_modesKept<=0 and upstream_coil_compression_thres>0, this threshold will be used to determine how many channels to keep
        /// the first N and first S will be used to compute number of channels to keep
        GADGET_PROPERTY(upstream_coil_compression_thres, double, "Threadhold for upstream coil compression", -1);
        GADGET_PROPERTY(upstream_coil_compression_num_modesKept, int, "Number of modes to keep for upstream coil compression", 0);

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        /// ------------------------------------------------------------------------------------
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // number of encoding spaces in the protocol
        size_t num_encoding_spaces_;

        // for every encoding space
        // calibration mode
        std::vector<Gadgetron::ismrmrdCALIBMODE> calib_mode_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // number of times the process function is called
        size_t process_called_times_;

        // store the KLT coefficients for N, S, SLC at every encoding space
        std::vector< std::vector< std::vector< std::vector< KLTType > > > > KLT_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // in verbose mode, more info is printed out
        bool verbose_;

        //// debug folder
        // std::string debug_folder_full_path_;

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        //// exporter
        // Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);
    };
}
