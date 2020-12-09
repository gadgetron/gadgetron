//
// Created by dchansen on 2/21/19.
//

#pragma once

#include <boost/range/adaptor/strided.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>
#include <range/v3/action.hpp>

#include "PureGadget.h"
#include "cmr_kspace_binning.h"
#include "mri_core_data.h"
#include "ImageArraySendMixin.h"

namespace Gadgetron {
    class PureCmrCartesianKSpaceBinningCineGadget : public Core::PureGadget<IsmrmrdImageArray, IsmrmrdReconData>, public ImageArraySendMixin<PureCmrCartesianKSpaceBinningCineGadget> {
    public:
        PureCmrCartesianKSpaceBinningCineGadget(const Core::Context& context, const Core::GadgetProperties& props);

        IsmrmrdImageArray process_function(IsmrmrdReconData args) const override;

        NODE_PROPERTY(verbose,bool,"Verbose",false);
        /// parameters for workflow
        NODE_PROPERTY(
            use_multiple_channel_recon, bool, "Whether to perform multi-channel recon in the raw data step", true);
        NODE_PROPERTY(use_nonlinear_binning_recon, bool, "Whether to non-linear recon in the binning step", true);
        NODE_PROPERTY(number_of_output_phases, size_t, "Number of output phases after binning", 30);

        NODE_PROPERTY(send_out_raw, bool, "Whether to set out raw images", false);
        NODE_PROPERTY(
            send_out_multiple_series_by_slice, bool, "Whether to set out binning images as multiple seires", false);

        /// parameters for raw image reconstruction
        NODE_PROPERTY(arrhythmia_rejector_factor, float,
            "If a heart beat RR is not in the range of [ (1-arrhythmiaRejectorFactor)*meanRR "
            "(1+arrhythmiaRejectorFactor)*meanRR], it will be rejected",
            0.25);

        NODE_PROPERTY(grappa_kSize_RO, int, "Raw data recon, kernel size RO", 5);
        NODE_PROPERTY(grappa_kSize_E1, int, "Raw data recon, kernel size E1", 4);
        NODE_PROPERTY(grappa_reg_lamda, double, "Raw data recon, kernel calibration regularization", 0.0005);

        NODE_PROPERTY(downstream_coil_compression_num_modesKept, size_t,
            "Number of modes kept for down stream coil compression in raw recon step", 0);
        NODE_PROPERTY(downstream_coil_compression_thres, double,
            "Threshold for determining number of kept modes in the down stream coil compression", 0.01);

        /// parameters for kspace binning
        NODE_PROPERTY(
            respiratory_navigator_moco_reg_strength, double, "Regularization strength of respiratory moco", 6.0);

        NODE_PROPERTY(respiratory_navigator_moco_iters, std::vector<unsigned int>,
            "Number of iterations for respiratory moco", (std::vector<unsigned int>{ 1, 32, 100, 100 }));

        NODE_PROPERTY(
            kspace_binning_interpolate_heart_beat_images, bool, "Whether to interpolate best heart beat images", true);
        NODE_PROPERTY(
            kspace_binning_navigator_acceptance_window, double, "Respiratory navigator acceptance window", 0.65);
        NODE_PROPERTY(kspace_binning_max_temporal_window, double,
            "Maximally allowed temporal window ratio for binned kspace", 2.0);
        NODE_PROPERTY(kspace_binning_minimal_cardiac_phase_width, double,
            "Allowed  minimal temporal window for binned kspace, in ms", 25.0);

        NODE_PROPERTY(kspace_binning_moco_reg_strength, double, "Regularization strength of binning moco", 12.0);
        NODE_PROPERTY(kspace_binning_moco_iters, std::vector<unsigned int>, "Number of iterations for binning moco",
            (std::vector<unsigned int>{ 24, 64, 100, 100, 100 }));

        /// parameters for recon on binned kspace
        NODE_PROPERTY(kspace_binning_kSize_RO, int, "Binned kspace recon, kernel size RO", 7);
        NODE_PROPERTY(kspace_binning_kSize_E1, int, "Binned kspace recon, kernel size E1", 7);
        NODE_PROPERTY(
            kspace_binning_reg_lamda, double, "Binned kspace recon, kernel calibration regularization", 0.005);
        /// linear recon step
        NODE_PROPERTY(kspace_binning_linear_iter_max, size_t,
            "Binned kspace recon, maximal number of iterations, linear recon", 90);
        NODE_PROPERTY(
            kspace_binning_linear_iter_thres, double, "Binned kspace recon, iteration threshold, linear recon", 0.0015);
        /// Non-linear recon step
        NODE_PROPERTY(kspace_binning_nonlinear_iter_max, size_t,
            "Binned kspace recon, maximal number of iterations, non-linear recon", 25);
        NODE_PROPERTY(kspace_binning_nonlinear_iter_thres, double,
            "Binned kspace recon, iteration threshold, non-linear recon", 0.004);
        NODE_PROPERTY(kspace_binning_nonlinear_data_fidelity_lamda, double,
            "Binned kspace recon, strength of data fidelity term, non-linear recon", 1.0);
        NODE_PROPERTY(kspace_binning_nonlinear_image_reg_lamda, double,
            "Binned kspace recon, strength of image term, non-linear recon", 0.00015);
        NODE_PROPERTY(kspace_binning_nonlinear_reg_N_weighting_ratio, double,
            "Binned kspace recon, regularization weighting ratio along the N dimension, non-linear recon", 10.0);
        NODE_PROPERTY(kspace_binning_nonlinear_reg_use_coil_sen_map, bool,
            "Binned kspace recon, whether to use coil map, non-linear recon", false);
        NODE_PROPERTY(kspace_binning_nonlinear_reg_with_approx_coeff, bool,
            "Binned kspace recon, whether to keep approximal coefficients, non-linear recon", true);
        NODE_PROPERTY(kspace_binning_nonlinear_reg_wav_name, std::string,
            "Binned kspace recon, wavelet name, non-linear recon. Possible values: db1, db2, db3, db4, db5", "db1");
        NODE_PROPERTY(time_tick, float, "Time tick in ms", 2.5f);

    private:
        CmrKSpaceBinning<float> create_binner() const;

        struct BinningResult {
            IsmrmrdImageArray image;
            hoNDArray<float> acquisition_time;
            hoNDArray<float> capture_time;
        };

        // acceleration factor for E1 and E2
        std::vector<double> acceFactorE1_;
        std::vector<double> acceFactorE2_;
        // kspace line offset in case of partial fourier
        std::vector<int> space_matrix_offset_E1_;
        std::vector<int> space_matrix_offset_E2_;

        // calibration mode
        std::vector<Gadgetron::ismrmrdCALIBMODE> calib_mode_;

        BinningResult perform_binning(IsmrmrdReconBit reconBit, size_t encoding) const;
        void set_image_header(const IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t enc) const;
    };

}
