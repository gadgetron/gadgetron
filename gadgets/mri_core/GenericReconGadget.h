/** \file   GenericReconGadget.h
    \brief  This serves an optional base class gadget for both 2DT and 3DT reconstruction, working on the IsmrmrdReconData.
            Some common functionalities are implemented here and can be reused in specific recon gadgets.
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
#include "hoNDFFT.h"

#include "mri_core_def.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"
#include "mri_core_coil_map_estimation.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconGadget : public Gadget1<IsmrmrdReconData>
    {
    public:
        GADGET_DECLARE(GenericReconGadget);

        typedef Gadget1<IsmrmrdReconData> BaseClass;

        GenericReconGadget();
        ~GenericReconGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// image series
        GADGET_PROPERTY(image_series, int, "Image series number", 0);

        /// coil map estimation method
        GADGET_PROPERTY_LIMITS(coil_map_algorithm, std::string, "Method for coil map estimation", "Inati",
            GadgetPropertyLimitsEnumeration, "Inati", "Inati_Iter");

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

        // for every encoding space

        // acceleration factor for E1 and E2
        std::vector<double> acceFactorE1_;
        std::vector<double> acceFactorE2_;

        // calibration mode
        std::vector<Gadgetron::ismrmrdCALIBMODE> calib_mode_;

        // encoding space limits
        std::vector<ISMRMRD::EncodingCounters> meas_max_idx_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------
        // number of times the process function is called
        size_t process_called_times_;

        /// buffers used during recon
        hoNDArray< std::complex<float> > complex_im_recon_buf_;
        hoNDArray< std::complex<float> > data_recon_buf_;

        // filter used for ref coil map 
        hoNDArray< std::complex<float> > filter_RO_ref_coi_map_;
        hoNDArray< std::complex<float> > filter_E1_ref_coi_map_;
        hoNDArray< std::complex<float> > filter_E2_ref_coi_map_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        //// debug folder
        //std::string debug_folder_full_path_;

        //// exporter
        //Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);

        // --------------------------------------------------
        // recon step functions
        // --------------------------------------------------

        // make the ref data for coil map estimation
        virtual void make_ref_coil_map(IsmrmrdDataBuffered& ref_, std::vector<size_t> recon_dims, hoNDArray< std::complex<float> >& ref_calib, hoNDArray< std::complex<float> >& ref_coil_map, size_t encoding);

        // estimate coil map
        virtual void perform_coil_map_estimation(const hoNDArray< std::complex<float> >& ref_coil_map, hoNDArray< std::complex<float> >& coil_map, size_t encoding);

        // compute image header
        virtual void compute_image_header(IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t encoding);

        // send out the recon results
        virtual int send_out_image_array(IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t encoding, int series_num, const std::string& data_role);

        // --------------------------------------------------
        // utility functions
        // --------------------------------------------------
        // compute image number
        virtual size_t compute_image_number(ISMRMRD::ImageHeader& imheader, size_t encoding = 0, size_t CHA = 1, size_t cha = 0, size_t E2 = 1);
    };
}
