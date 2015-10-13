/** \file   GenericReconFieldOfViewAdjustmentGadget.h
    \brief  This is the class gadget for both 2DT and 3DT reconstruction, working on the IsmrmrdImageArray.
            This gadget will adjust the image field-of-view and/or image size to the protocol prescribed values.

            This class is a part of general cartesian recon chain.

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

    class EXPORTGADGETSMRICORE GenericReconFieldOfViewAdjustmentGadget : public Gadget1<IsmrmrdImageArray>
    {
    public:
        GADGET_DECLARE(GenericReconFieldOfViewAdjustmentGadget);

        typedef Gadget1<IsmrmrdImageArray> BaseClass;

        GenericReconFieldOfViewAdjustmentGadget();
        ~GenericReconFieldOfViewAdjustmentGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

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
        // recon size
        std::vector< std::vector<size_t> > recon_size_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // kspace filter
        hoNDArray< std::complex<float> > filter_RO_;
        hoNDArray< std::complex<float> > filter_E1_;
        hoNDArray< std::complex<float> > filter_E2_;

        // kspace buffer
        hoNDArray< std::complex<float> > kspace_buf_;

        // results of filtering
        hoNDArray< std::complex<float> > res_;

        // number of times the process function is called
        size_t process_called_times_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // debug folder
        // std::string debug_folder_full_path_;

        // clock for timing
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

        // adjust FOV
        int adjust_FOV(IsmrmrdImageArray& data);

        // perform fft or ifft
        void perform_fft(size_t E2, const hoNDArray< std::complex<float> >& input, hoNDArray< std::complex<float> >& output);
        void perform_ifft(size_t E2, const hoNDArray< std::complex<float> >& input, hoNDArray< std::complex<float> >& output);
    };
}
