/** \file   GenericReconGadget.h
    \brief  This serves an optional base class gadget for both 2DT and 3DT reconstruction, working on the IsmrmrdReconData.
            Some common functionalities are implemented here and can be reused in specific recon gadgets.
    \author Hui Xue
*/

#pragma once

#include "GenericReconBase.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDFFT.h"

#include "mri_core_coil_map_estimation.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconGadget : public GenericReconDataBase
    {
    public:
        GADGET_DECLARE(GenericReconGadget);

        typedef GenericReconDataBase BaseClass;

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

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

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

        /// buffers used during recon
        hoNDArray< std::complex<float> > complex_im_recon_buf_;
        hoNDArray< std::complex<float> > data_recon_buf_;

        // filter used for ref coil map 
        hoNDArray< std::complex<float> > filter_RO_ref_coi_map_;
        hoNDArray< std::complex<float> > filter_E1_ref_coi_map_;
        hoNDArray< std::complex<float> > filter_E2_ref_coi_map_;

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

        // prepare header for sending out
        virtual int prep_image_header_send_out(IsmrmrdImageArray& res, size_t n, size_t s, size_t slc, size_t encoding, int series_num, const std::string& data_role);

        // send out the recon results
        virtual int send_out_image_array(IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t encoding, int series_num, const std::string& data_role);

        // compute snr scaling factor from effective acceleration rate and sampling region
        void compute_snr_scaling_factor(IsmrmrdReconBit& recon_bit, float& effective_acce_factor, float& snr_scaling_ratio);

        // --------------------------------------------------
        // utility functions
        // --------------------------------------------------
        // compute image number
        virtual size_t compute_image_number(ISMRMRD::ImageHeader& imheader, size_t encoding = 0, size_t CHA = 1, size_t cha = 0, size_t E2 = 1);
    };
}
