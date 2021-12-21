/** \file   CmrCartesianKSpaceBinningCineGadget.h
    \brief  This is the class gadget for 2DT cartesian KSpace binning reconstruction, working on the IsmrmrdReconData.
            The temporal dimension is dimension N and motion sharing dimension is S.

            Ref: 
            [1] Hui Xue, Peter Kellman, Gina LaRocca, Andrew E Arai and Michael S Hansen. 
                High spatial and temporal resolution retrospective cine cardiovascular magnetic resonance from shortened free breathing real-time acquisitions.
                Journal of Cardiovascular Magnetic Resonance 2013; 15:102.

            [2] Peter Kellman, Christophe Chefd’hotel, Christine H. Lorenz, Christine Mancini, Andrew E. Arai, Elliot R. McVeigh. 
                High spatial and temporal resolution cardiac cine MRI from retrospective reconstruction of data acquired in real time using motion correction and resorting. 
                Magn Reson Med. 2009; 62:1557–64.

            [3] Michael S. Hansen, Thomas S. Sørensen, Andrew E. Arai, and Peter Kellman. 
                Retrospective reconstruction of high temporal resolution cine images from real-time MRI using iterative motion correction. 
                Magn Reson Med. 2012; 68:741–50.

    \author Hui Xue
*/

#pragma once

#include "gadgetron_cmr_export.h"
#include "GenericReconGadget.h"
#include "cmr_kspace_binning.h"

namespace Gadgetron {

    class EXPORTGADGETSCMR CmrCartesianKSpaceBinningCineGadget : public GenericReconGadget
    {
    public:
        GADGET_DECLARE(CmrCartesianKSpaceBinningCineGadget);

        typedef GenericReconGadget BaseClass;

        CmrCartesianKSpaceBinningCineGadget();
        ~CmrCartesianKSpaceBinningCineGadget();

        /// parameters for workflow
        GADGET_PROPERTY(use_multiple_channel_recon, bool, "Whether to perform multi-channel recon in the raw data step", true);
        GADGET_PROPERTY(use_nonlinear_binning_recon, bool, "Whether to non-linear recon in the binning step", true);
        GADGET_PROPERTY(number_of_output_phases, int, "Number of output phases after binning", 30);

        GADGET_PROPERTY(send_out_raw, bool, "Whether to set out raw images", false);
        GADGET_PROPERTY(send_out_multiple_series_by_slice, bool, "Whether to set out binning images as multiple series", false);

        /// parameters for raw image reconstruction
        GADGET_PROPERTY(arrhythmia_rejector_factor, float, "If a heart beat RR is not in the range of [ (1-arrhythmiaRejectorFactor)*meanRR (1+arrhythmiaRejectorFactor)*meanRR], it will be rejected", 0.25);

        GADGET_PROPERTY(grappa_kSize_RO, int, "Raw data recon, kernel size RO", 5);
        GADGET_PROPERTY(grappa_kSize_E1, int, "Raw data recon, kernel size E1", 4);
        GADGET_PROPERTY(grappa_reg_lamda, double, "Raw data recon, kernel calibration regularization", 0.0005);

        GADGET_PROPERTY(downstream_coil_compression_num_modesKept, size_t, "Number of modes kept for down stream coil compression in raw recon step", 0);
        GADGET_PROPERTY(downstream_coil_compression_thres, double, "Threshold for determining number of kept modes in the down stream coil compression", 0.01);

        /// parameters for kspace binning
        GADGET_PROPERTY(respiratory_navigator_moco_reg_strength, double, "Regularization strength of respiratory moco", 6.0);
        GadgetProperty<std::vector<unsigned int>, GadgetPropertyLimitsNoLimits<std::vector<unsigned int> > > respiratory_navigator_moco_iters{"respiratory_navigator_moco_iters", 
                                                                                                            "unsigned int", 
                                                                                                            "Number of iterations for respiratory moco", 
                                                                                                            this, 
                                                                                                            {1, 32, 100, 100}, 
                                                                                                            GadgetPropertyLimitsNoLimits<std::vector<unsigned int> >()};

        GADGET_PROPERTY(kspace_binning_interpolate_heart_beat_images, bool, "Whether to interpolate best heart beat images", true);
        GADGET_PROPERTY(kspace_binning_navigator_acceptance_window, double, "Respiratory navigator acceptance window", 0.65);
        GADGET_PROPERTY(kspace_binning_max_temporal_window, double, "Maximally allowed temporal window ratio for binned kspace", 2.0);
        GADGET_PROPERTY(kspace_binning_minimal_cardiac_phase_width, double, "Allowed  minimal temporal window for binned kspace, in ms", 25.0);

        GADGET_PROPERTY(kspace_binning_moco_reg_strength, double, "Regularization strength of binning moco", 12.0);
        GadgetProperty<std::vector<unsigned int>, GadgetPropertyLimitsNoLimits<std::vector<unsigned int> > > kspace_binning_moco_iters{"kspace_binning_moco_iters", 
                                                                                                            "unsigned int", 
                                                                                                            "Number of iterations for binning moco", 
                                                                                                            this, 
                                                                                                            {24, 64, 100, 100, 100}, 
                                                                                                            GadgetPropertyLimitsNoLimits<std::vector<unsigned int> >()};

        /// parameters for recon on binned kspace
        GADGET_PROPERTY(kspace_binning_kSize_RO, int, "Binned kspace recon, kernel size RO", 7);
        GADGET_PROPERTY(kspace_binning_kSize_E1, int, "Binned kspace recon, kernel size E1", 7);
        GADGET_PROPERTY(kspace_binning_reg_lamda, double, "Binned kspace recon, kernel calibration regularization", 0.005);
        /// linear recon step
        GADGET_PROPERTY(kspace_binning_linear_iter_max, size_t, "Binned kspace recon, maximal number of iterations, linear recon", 90);
        GADGET_PROPERTY(kspace_binning_linear_iter_thres, double, "Binned kspace recon, iteration threshold, linear recon", 0.0015);
        /// Non-linear recon step
        GADGET_PROPERTY(kspace_binning_nonlinear_iter_max, size_t, "Binned kspace recon, maximal number of iterations, non-linear recon", 25);
        GADGET_PROPERTY(kspace_binning_nonlinear_iter_thres, double, "Binned kspace recon, iteration threshold, non-linear recon", 0.004);
        GADGET_PROPERTY(kspace_binning_nonlinear_data_fidelity_lamda, double, "Binned kspace recon, strength of data fidelity term, non-linear recon", 1.0);
        GADGET_PROPERTY(kspace_binning_nonlinear_image_reg_lamda, double, "Binned kspace recon, strength of image term, non-linear recon", 0.00015);
        GADGET_PROPERTY(kspace_binning_nonlinear_reg_N_weighting_ratio, double, "Binned kspace recon, regularization weighting ratio along the N dimension, non-linear recon", 10.0);
        GADGET_PROPERTY(kspace_binning_nonlinear_reg_use_coil_sen_map, bool, "Binned kspace recon, whether to use coil map, non-linear recon", false);
        GADGET_PROPERTY(kspace_binning_nonlinear_reg_with_approx_coeff, bool, "Binned kspace recon, whether to keep approximal coefficients, non-linear recon", true);
        GADGET_PROPERTY_LIMITS(kspace_binning_nonlinear_reg_wav_name, std::string, "Binned kspace recon, wavelet name, non-linear recon", "db1", 
            GadgetPropertyLimitsEnumeration, "db1", "db2", "db3", "db4", "db5");

        // for debug
        GadgetProperty<std::vector<unsigned int>, GadgetPropertyLimitsNoLimits<std::vector<unsigned int> > > kspace_binning_processed_slices{"kspace_binning_processed_slices", 
                                                                                                            "unsigned int", 
                                                                                                            "If set, these slices will be processed", 
                                                                                                            this, 
                                                                                                            {}, 
                                                                                                            GadgetPropertyLimitsNoLimits<std::vector<unsigned int> >()};

    protected:

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // binning object
        Gadgetron::CmrKSpaceBinning<float> binning_reconer_;

        // the raw recon results
        // [RO E1 E2 1 N S SLC]
        IsmrmrdImageArray res_raw_;
        // acqusition time  and trigger time in ms
        // [N S SLC]
        hoNDArray< float > acq_time_raw_;
        hoNDArray< float > cpt_time_raw_;

        // the binning recon results
        // [RO E1 E2 1 binned_N S SLC]
        IsmrmrdImageArray res_binning_;
        // acqusition time  and trigger time in ms for binned images
        // [binned_N S SLC]
        hoNDArray< float > acq_time_binning_;
        hoNDArray< float > cpt_time_binning_;

        // if true, every slice will be sent out as separate series
        bool send_out_multiple_series_by_slice_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);

        // --------------------------------------------------
        // recon step functions
        // --------------------------------------------------
        virtual void perform_binning(IsmrmrdReconBit& recon_bit, size_t encoding);

        // create binning image header
        void create_binning_image_headers_from_raw();

        // set the time stamps
        void set_time_stamps(IsmrmrdImageArray& res, hoNDArray< float >& acq_time, hoNDArray< float >& cpt_time);

        // --------------------------------------------------
        // overload functions
        // --------------------------------------------------
        // send out the recon results
        virtual int prep_image_header_send_out(IsmrmrdImageArray& res, size_t n, size_t s, size_t slc, size_t encoding, int series_num, const std::string& data_role);
    };
}
