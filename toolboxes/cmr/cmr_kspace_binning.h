/** \file   cmr_kspace_binning.h
    \brief  Implement kspace binning recon for 2D acquisition
            The input has dimension [RO E1 CHA N S SLC]
            Temporal dimension is N and S is motion sharing dimension, e.g. contrast or set
    \author Hui Xue
*/

#pragma once

#include "cmr_export.h"

#include "GadgetronTimer.h"

#ifdef min
#undef min
#endif // min

#include <algorithm>
#include "hoMatrix.h"

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"

#include "mri_core_def.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"

#include "ImageIOAnalyze.h"
#include "hoImageRegContainer2DRegistration.h"

namespace Gadgetron { 

    template <typename T>
    class EXPORTCMR KSpaceBinningObj
    {
    public:

        KSpaceBinningObj();
        virtual ~KSpaceBinningObj();

        // ------------------------------------
        /// recon outputs
        // ------------------------------------
        /// binning reconstructed images, headers and meta attributes
        /// [RO E1 1 output_N_ S]
        hoNDArray< std::complex<T> > complex_image_binning_;

        /// raw reconstructed images, headers and meta attributes
        /// [RO E1 dstCHA N S]
        hoNDArray< std::complex<T> > full_kspace_raw_;
        /// complex images, [RO E1 1 N S]
        hoNDArray< std::complex<T> > complex_image_raw_;
        /// coil map, [RO E1 dstCHA S]
        hoNDArray< std::complex<T> > coil_map_raw_;

        /// mean RR interval in ms
        float mean_RR_;

        // ------------------------------------
        /// recon inputs
        // ------------------------------------
        /// number of required binning output along N
        size_t output_N_;

        /// acceleration factor
        float accel_factor_E1_;

        /// whether random sampling is used
        bool random_sampling_;

        /// input kspace data and headers
        /// data: [RO E1 CHA N S]
        /// only 2D binning is supported, so E2==1
        /// N is the temporal dimension
        /// S is the motion sharing dimension, e.g. contrast or set
        /// incoming image should already be in eigeh channel, with the high energy channel first
        hoNDArray< std::complex<T> > data_;
        SamplingDescription sampling_;

        /// [E1 N S]
        hoNDArray< ISMRMRD::AcquisitionHeader > headers_;

        // ------------------------------------
        /// buffer for recon
        // ------------------------------------
        /// time stamp in ms for kspace lines [E1 N S]
        hoNDArray<float> time_stamp_;
        /// cardiac phase time in ms for kspace lines, [E1 N S]
        hoNDArray<float> cpt_time_stamp_;
        /// cardiac phase time ratio (CPT) [0 1] for kspace lines, [E1 N S]
        hoNDArray<float> cpt_time_ratio_;

        /// time stamp for every raw images, [N S]
        hoNDArray<float> phs_time_stamp_;
        /// cardiac phase time in ms for every raw images, [N S]
        hoNDArray<float> phs_cpt_time_stamp_;
        /// cardiac phase time ratio for every raw images, [N S]
        hoNDArray<float> phs_cpt_time_ratio_;

        // buffer to store index of heart beat for every kspace line; the first heart beat has index 0
        hoNDArray<int> ind_heart_beat_;
        // arrays to store the starting and ending [e1 n] indexes for every heart beat
        // for every S, one std::vector is created to store starting/ending of heart beats
        std::vector< std::vector< std::pair<size_t, size_t> > > starting_heart_beat_;
        std::vector< std::vector< std::pair<size_t, size_t> > > ending_heart_beat_;

        // output cardiac phase time stamp in ms
        std::vector<float> desired_cardiac_phases_;
        // output cardiac phase ratio
        std::vector<float> desired_cpt_;

        /// respiratory navigator signal for every kspace line, [E1 N S]
        hoNDArray<float> navigator_;

        /// binned kspace
        /// [RO E1 1ordstCHA output_N_ S]
        hoNDArray< std::complex<T> > kspace_binning_;
        hoNDArray< std::complex<T> > kspace_binning_wider_;
        /// image domain averaged kspace, [RO E1 1ordstCHA output_N_ S]
        hoNDArray< std::complex<T> > kspace_binning_image_domain_average_;
        /// the number of lines falling into the binned kspace
        /// [E1 output_N_ S]
        hoNDArray< float > kspace_binning_hit_count_;
    };

    template<typename T> 
    class CmrKSpaceBinning
    {
    public:

        typedef CmrKSpaceBinning<T> Self;
        typedef hoNDArray< std::complex<T> > ArrayType;

        // for every heart beat, record its staring and ending [e1 n] indexes
        typedef std::vector< std::pair<size_t, size_t> > HeartBeatIndexType;

        // store the respiratory navigator detected ROI
        typedef std::vector< std::pair< hoNDPoint<float, 2>, hoNDPoint<float, 2> > > NavigatorRoiType;

        // image type
        typedef Gadgetron::hoNDImage<T, 2> ImageType;
        typedef Gadgetron::hoNDImageContainer2D<ImageType> ImageContinerType;

        // image type
        typedef Gadgetron::hoNDImage< std::complex<T>, 2> ComplexImageType;
        typedef Gadgetron::hoNDImageContainer2D<ComplexImageType> ComplexImageContinerType;

        // motion correction
        typedef Gadgetron::hoImageRegContainer2DRegistration<ImageType, ImageType, double> RegContainer2DType;

        // deformation field
        typedef hoNDImageContainer2D< hoNDImage<double, 2> > DeformationFieldContinerType;

        CmrKSpaceBinning();
        virtual ~CmrKSpaceBinning();

        // ======================================================================================
        // perform the binning reconstruction
        // ======================================================================================
        virtual void process_binning_recon();

        // ------------------------------------
        /// binning object, storing the kspace data and results
        // ------------------------------------
        KSpaceBinningObj<T> binning_obj_;

        // ======================================================================================
        /// parameter for overall workflow
        // ======================================================================================

        // whether true, raw data recon step will produce multiple channels
        bool use_multiple_channel_recon_;

        // if true, and use_multiple_channel_recon_==true, then
        // the binning recon step will perform parallel imaging
        bool use_paralell_imaging_binning_recon_;

        // if true, and use_multiple_channel_recon_==ture and use_paralell_imaging_binning_recon_==true
        // the binning recon step will perform nonlinear reconstruction
        bool use_nonlinear_binning_recon_;

        // if true, respiratory navigator will be estimated
        // if false, user can supply external navigator by filling in binning_obj_.navigator_
        bool estimate_respiratory_navigator_;

        // ======================================================================================
        /// parameter for respiratory navigator estimation
        // ======================================================================================

        // regularization strength of respiratory navigator moco
        T respiratory_navigator_moco_reg_strength_;

        // number of iterations
        std::vector<unsigned int> respiratory_navigator_moco_iters_;

        // respiratory detection patch size and step size
        size_t respiratory_navigator_patch_size_RO_;
        size_t respiratory_navigator_patch_size_E1_;

        size_t respiratory_navigator_patch_step_size_RO_;
        size_t respiratory_navigator_patch_step_size_E1_;

        // ======================================================================================
        /// parameter for raw image reconstruction
        // ======================================================================================

        // time tick for every time stamp unit, in ms (e.g. 2.5ms)
        float time_tick_;

        // index of cardiac trigger time in the physilogical time stamp array
        size_t trigger_time_index_;

        // if a heart beat RR is not in the range of [ (1-arrhythmiaRejectorFactor)*meanRR (1+arrhythmiaRejectorFactor)*meanRR], it will be rejected
        float arrhythmia_rejector_factor_;

        // ----------------------
        // grappa raw data reconstruction
        // ----------------------
        int grappa_kSize_RO_;
        int grappa_kSize_E1_;
        double grappa_reg_lamda_;

        // down stream coil compression
        // if downstream_coil_compression_num_modesKept > 0, this number of channels will be used as the dst channels
        // if downstream_coil_compression_num_modesKept==0 and downstream_coil_compression_thres>0, 
        // the number of dst channels will be determined  by this threshold
        size_t downstream_coil_compression_num_modesKept_;
        double downstream_coil_compression_thres_;

        // ----------------------
        // spirit raw data reconstruction, if random sampling was used
        // ----------------------
        // to be implemented

        // ======================================================================================
        /// parameter for cardiac binning
        // ======================================================================================

        // whether to interpolate heart beat images; if not, interpolate the deformation fields
        bool kspace_binning_interpolate_heart_beat_images_;

        // acceptance window for respiratory navigator
        float kspace_binning_navigator_acceptance_window_;

        // regularization strength of binning moco
        T kspace_binning_moco_reg_strength_;

        // number of iterations for binning moco
        std::vector<unsigned int> kspace_binning_moco_iters_;

        // max allowed temporal ratio window
        float kspace_binning_max_temporal_window_;

        // minimal allowed temporal window in ms
        float kspace_binning_minimal_cardiac_phase_width_;

        // ======================================================================================
        /// parameter for recon after binning
        // ======================================================================================

        // kernel size and calibration regularization threshold on the binned kspace
        int kspace_binning_kSize_RO_;
        int kspace_binning_kSize_E1_;
        double kspace_binning_reg_lamda_;

        // maximal number of iterations for linear recon on binned kspace
        size_t kspace_binning_linear_iter_max_;
        // threshold to stop iterations for linear recon on binned kspace
        double kspace_binning_linear_iter_thres_;

        // maximal number of iterations for non-linear recon on binned kspace
        size_t kspace_binning_nonlinear_iter_max_;
        // threshold to stop the non-linear iterations
        double kspace_binning_nonlinear_iter_thres_;
        // data fidelity term; if 0, null-space constraint is applied
        double kspace_binning_nonlinear_data_fidelity_lamda_;
        // image regularization term
        double kspace_binning_nonlinear_image_reg_lamda_;
        // regularization weighting ratio along the N dimension
        double kspace_binning_nonlinear_reg_N_weighting_ratio_;
        // whether to use coil map in the image regularization term
        bool kspace_binning_nonlinear_reg_use_coil_sen_map_;
        // whether to keep approximation coefficients during regularization
        bool kspace_binning_nonlinear_reg_with_approx_coeff_;
        // wavelet type for nonlinear reg
        std::string kspace_binning_nonlinear_reg_wav_name_;

        // ======================================================================================
        /// parameter for debugging
        // ======================================================================================
        std::string suffix_;
        bool verbose_;
        std::string debug_folder_;
        bool perform_timing_;

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        // exporter
        Gadgetron::ImageIOAnalyze gt_exporter_;

    protected:

        // ======================================================================================
        // perform every steps
        // ======================================================================================

        // perform raw data recon on binning_obj_.data_
        // fill the binning_obj_.full_kspace_raw_, binning_obj_.complex_image_raw_, binning_obj_.coil_map_raw_
        virtual void perform_raw_data_recon();

        /// this function parses all recorded time stamps and computes cpt time stamps
        /// first, all time stamps are converted into the unit of ms
        /// second, a linear fit is performed to reduce then quantization error
        /// third, for all missing lines, its time stamps are estimated
        /// fourth, the cpt time stamps are analyzed to find out every heart beats
        /// fifth, the mean RR and desired cardiac phases are computed
        virtual void estimate_time_stamps();

        /// estimate the respiratory navigator from the images
        virtual void estimate_respiratory_navigator();

        /// find the best heart beat
        /// the strategy is to find the heart beat with minimal intra-heart-beat navigator signal variance
        /// bestHB for every S
        virtual void find_best_heart_beat(std::vector<size_t>& bestHB);

        /// reject heart beat with irregular rhythm
        virtual void reject_irregular_heart_beat();

        /// compute kspace after binning
        virtual void compute_kspace_binning(const std::vector<size_t>& bestHB, std::vector<size_t>& slices_not_processing);

        /// perform recon on the binned kspace
        virtual void perform_recon_binned_kspace(const std::vector<size_t>& slices_not_processing);

        // ======================================================================================
        // implementation functions
        // ======================================================================================
        /// if the alternativing acqusition is used, detect and flip the time stamps
        void detect_and_flip_alternating_order(const hoNDArray<float>& time_stamp, const hoNDArray<float>& cpt_time_stamp, hoNDArray<float>& cpt_time_stamp_flipped, std::vector<bool>& ascending);

        /// from the input time_stamp and cpt_time_stamp, compute times for missing lines, find completed heart beat and compute times for every phase (n)
        /// time_stamp : [E1 N]
        /// cpt_time_stamp : [E1 N]
        void process_time_stamps(hoNDArray<float>& time_stamp, hoNDArray<float>& cpt_time_stamp, 
                                hoNDArray<float>& cpt_time_ratio, hoNDArray<float>& phs_time_stamp, 
                                hoNDArray<float>& phs_cpt_time_stamp, hoNDArray<float>& phs_cpt_time_ratio, 
                                hoNDArray<int>& indHeartBeat, HeartBeatIndexType& startingHB, 
                                HeartBeatIndexType& endingHB, float& meanRRInterval);

        /// compute RR interval for a heart beat
        void compute_RRInterval(size_t s, size_t HB, float& RRInterval);

        /// compute navigator metrics for a heart beat
        void compute_metrics_navigator_heart_beat(size_t s, size_t HB, float& mean_resp, float& var_resp, float& min_resp, float& max_resp);

        /// interpolate the best heart beat images to get images at the desired cardiace phase ration
        void interpolate_best_HB_images(const std::vector<float>& cpt_time_ratio_bestHB, const hoNDArray<T>& mag_bestHB, const std::vector<float>& desired_cpt, hoNDArray<T>& mag_bestHB_at_desired_cpt);

        /// select images using the navigator and cpt time ratio
        /// selected: store the selected image n indexes
        void select_images_with_navigator(float desired_cpt_time_ratio, float accepted_nav[2], size_t n, size_t s, std::vector<size_t>& selected);

        /// perform moco between best heart beat images and selected images
        void perform_moco_selected_images_with_best_heart_beat(const std::vector < std::vector<size_t> >& selected_images, const hoNDArray<T>& mag_s, const hoNDArray<T>& mag_bestHB_at_desired_cpt, DeformationFieldContinerType& dx, DeformationFieldContinerType& dy);

        /// warp complex images using the moco between best heart beat images and selected images
        void perform_moco_warp_on_selected_images(const std::vector < std::vector<size_t> >& selected_images, const ArrayType& complex_image, DeformationFieldContinerType deform[2], std::vector<ArrayType>& warpped_complex_images);

        /// fill the binned kspace
        void fill_binned_kspace(size_t s, size_t dst_n, const std::vector<size_t>& selected_images, const ArrayType& warpped_kspace, ArrayType& kspace_filled, hoNDArray<float>& hit_count);

        /// perform linear recon on binning
        void perform_linear_recon_on_kspace_binning(const ArrayType& underSampledKspace, const ArrayType& kspaceInitial, const ArrayType& coilMap, ArrayType& resKSpace, ArrayType& resIm, ArrayType& kernel, ArrayType& kernelIm);
        /// perform nonlinear recon for binning
        void perform_non_linear_recon_on_kspace_binning(const ArrayType& underSampledKspace, const ArrayType& kspaceLinear, const ArrayType& coilMap, const ArrayType& kernel, const ArrayType& kernelIm, ArrayType& resKSpace, ArrayType& resIm);
    };
}
