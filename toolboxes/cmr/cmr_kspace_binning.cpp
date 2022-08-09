/** \file   cmr_kspace_binning.cpp
    \brief  Implement kspace binning recon for 2D acquisition
    \author Hui Xue
*/

#include "cmr_kspace_binning.h"
#include "log.h"
#include "hoNDFFT.h"
#include "mri_core_grappa.h"
#include "mri_core_spirit.h"
#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_coil_map_estimation.h"

#include "cmr_time_stamp.h"
#include "cmr_motion_correction.h"
#include "cmr_spirit_recon.h"
#include <boost/math/special_functions/sign.hpp>

namespace Gadgetron { 

template <typename T> 
KSpaceBinningObj<T>::KSpaceBinningObj()
{
}

template <typename T> 
KSpaceBinningObj<T>::~KSpaceBinningObj()
{
}

// --------------------------------------------------------------------------

template <typename T> 
void simple_line_fit(const std::vector<T>& x, const std::vector<T>& y, T& a, T& b)
{
    try
    {
        size_t num = x.size();

        if(num<2)
        {
            a = 0;
            b = 0;
            return;
        }

        T sx(0), sy(0);

        size_t n;
        for (n=0; n<num; n++)
        {
            sx += x[n];
            sy += y[n];
        }

        T mx = sx / (T)(num);
        T syy = 0;
        b = 0;
        for (n=0; n<num; n++)
        {
            T v = (x[n] - mx);
            syy += v*v;
            b += v*y[n];
        }

        syy = (std::abs(syy) > 0 ? syy : boost::math::sign(syy)*FLT_EPSILON);
        b /= syy;
        a = (sy - sx*b) / (T)(num);
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in simple_line_fit ... ");
    }
}

template <typename T> 
CmrKSpaceBinning<T>::CmrKSpaceBinning()
{
    time_tick_ = 2.5;
    trigger_time_index_ = 0;

    arrhythmia_rejector_factor_ = 0.25;

    grappa_kSize_RO_ = 5;
    grappa_kSize_E1_ = 4;
    grappa_reg_lamda_ = 0.0005;

    downstream_coil_compression_num_modesKept_ = 0;
    downstream_coil_compression_thres_ = 0.005;

    respiratory_navigator_moco_reg_strength_ = 6.0;
    respiratory_navigator_moco_iters_.resize(4);
    respiratory_navigator_moco_iters_[0] = 1;
    respiratory_navigator_moco_iters_[1] = 32;
    respiratory_navigator_moco_iters_[2] = 100;
    respiratory_navigator_moco_iters_[3] = 100;

    respiratory_navigator_patch_size_RO_ = 32;
    respiratory_navigator_patch_size_E1_ = 32;

    respiratory_navigator_patch_step_size_RO_ = 10;
    respiratory_navigator_patch_step_size_E1_ = 10;

    kspace_binning_interpolate_heart_beat_images_ = true;

    kspace_binning_navigator_acceptance_window_ = 0.65;

    kspace_binning_moco_reg_strength_ = 12.0;

    kspace_binning_moco_iters_.resize(5, 0);
    kspace_binning_moco_iters_[0] = 24;
    kspace_binning_moco_iters_[1] = 64;
    kspace_binning_moco_iters_[2] = 100;
    kspace_binning_moco_iters_[3] = 100;
    kspace_binning_moco_iters_[4] = 100;

    kspace_binning_kSize_RO_ = 7;
    kspace_binning_kSize_E1_ = 7;
    kspace_binning_reg_lamda_ = 0.005;

    kspace_binning_linear_iter_max_ = 90;
    kspace_binning_linear_iter_thres_ = 0.0015;

    kspace_binning_nonlinear_iter_max_ = 15;
    kspace_binning_nonlinear_iter_thres_ = 0.002;
    kspace_binning_nonlinear_data_fidelity_lamda_ = 1.0;
    kspace_binning_nonlinear_image_reg_lamda_ = 0.00015;
    kspace_binning_nonlinear_reg_N_weighting_ratio_ = 10;
    kspace_binning_nonlinear_reg_use_coil_sen_map_ = false;
    kspace_binning_nonlinear_reg_with_approx_coeff_ = false;
    kspace_binning_nonlinear_reg_wav_name_ = "db1";

    verbose_ = false;
    perform_timing_ = false;

    gt_timer_.set_timing_in_destruction(false);
    gt_timer_local_.set_timing_in_destruction(false);
}

template <typename T> 
CmrKSpaceBinning<T>::~CmrKSpaceBinning()
{
}

template <typename T> 
void CmrKSpaceBinning<T>::process_binning_recon()
{
    try
    {
        hoNDArray< std::complex<T> >& kspace = binning_obj_.data_;
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t CHA = kspace.get_size(2);
        size_t N = kspace.get_size(3);
        size_t S = kspace.get_size(4);

        GDEBUG_STREAM("Processing suffix is " << suffix_)

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kspace, debug_folder_ + "binning_kspace" + suffix_);

        // -----------------------------------------------------
        // perform the raw data recon
        // -----------------------------------------------------
        if ( this->perform_timing_ ) { gt_timer_.start("raw data recon ... "); }

        this->perform_raw_data_recon();

        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(binning_obj_.full_kspace_raw_, debug_folder_ + "full_kspace_raw" + suffix_);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(binning_obj_.complex_image_raw_, debug_folder_ + "complex_image_raw" + suffix_);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(binning_obj_.coil_map_raw_, debug_folder_ + "coil_map_raw" + suffix_);

        // -----------------------------------------------------
        // estimate acquisition time stamp and cardiac phase time ratio for all acuqired kspace lines
        // count the heart beats and when it starts
        // -----------------------------------------------------
        binning_obj_.cpt_time_stamp_.create(E1, N, S);
        binning_obj_.ind_heart_beat_.create(E1, N, S);

        if ( this->perform_timing_ ) { gt_timer_.start("estimate time stamps ... "); }

        this->estimate_time_stamps();

        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.time_stamp_, debug_folder_ + "time_stamp" + suffix_);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.cpt_time_stamp_, debug_folder_ + "cpt_time_stamp" + suffix_);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.phs_time_stamp_, debug_folder_ + "phs_time_stamp" + suffix_);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.phs_cpt_time_stamp_, debug_folder_ + "phs_cpt_time_stamp" + suffix_);

        // -----------------------------------------------------
        // estimate respiratory nativagor signal
        // -----------------------------------------------------
        if ( this->perform_timing_ ) { gt_timer_.start("estimate respiratory navigator ... "); }

        if(estimate_respiratory_navigator_)
        {
            this->estimate_respiratory_navigator();
        }
        else
        {
            // if external navigator was not supplied, use all data
            if(binning_obj_.navigator_.get_size(0)!=E1 || binning_obj_.navigator_.get_size(1)!=N || binning_obj_.navigator_.get_size(2)!=S)
            {
                binning_obj_.navigator_.create(E1, N, S);
                Gadgetron::fill(binning_obj_.navigator_, 1.0f);
            }
        }
        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.navigator_, debug_folder_ + "respiratory_navigator" + suffix_);

        // -----------------------------------------------------
        // find the best heart beat from the respiratory navigator
        // -----------------------------------------------------
        std::vector<size_t> bestHB;
        this->find_best_heart_beat(bestHB);

        // -----------------------------------------------------
        // reject irregular heart beats
        // -----------------------------------------------------
        this->reject_irregular_heart_beat();

        // -----------------------------------------------------
        // release some memory to reduce peak RAM usage
        // -----------------------------------------------------
        if(binning_obj_.data_.delete_data_on_destruct()) binning_obj_.data_.clear();

        // -----------------------------------------------------
        // perform kspace binning
        // -----------------------------------------------------
        // all time stamps and raw full kspace is filled now
        // binning can be performed
        std::vector<size_t> slices_not_processing;
        this->compute_kspace_binning(bestHB, slices_not_processing);
        if(binning_obj_.full_kspace_raw_.delete_data_on_destruct()) binning_obj_.full_kspace_raw_.clear();

        // -----------------------------------------------------
        // perform recon on the binned kspace 
        // -----------------------------------------------------
        this->perform_recon_binned_kspace(slices_not_processing);
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::process_binning_recon() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::perform_raw_data_recon()
{
    try
    {
        // data: [RO E1 CHA N S]
        if(this->binning_obj_.random_sampling_)
        {
            GADGET_THROW("To be implemented, random sampling binning ... ");
        }

        ArrayType& data = this->binning_obj_.data_;
        // if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(data, debug_folder_ + "raw_data_recon_data" + suffix_);

        ArrayType& full_kspace = this->binning_obj_.full_kspace_raw_;
        ArrayType& complex_image = this->binning_obj_.complex_image_raw_;
        ArrayType& coil_map = this->binning_obj_.coil_map_raw_;

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t CHA = data.get_size(2);
        size_t N = data.get_size(3);
        size_t S = data.get_size(4);

        // ------------------------------------------------
        // average across N and S to compute ref
        // ------------------------------------------------
        bool count_sampling_freq = true;
        bool average_all_ref_N = true;
        bool average_all_ref_S = true;

        ArrayType data_for_ref(RO, E1, 1, CHA, N, S, 1, data.begin()), ref;

        if ( this->perform_timing_ ) { gt_timer_.start("--> perform_raw_data_recon, prepare ref"); }

        Gadgetron::compute_averaged_data_N_S(data_for_ref, average_all_ref_N, average_all_ref_S, count_sampling_freq, ref);

        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(ref, debug_folder_ + "raw_data_recon_ref" + suffix_);

        // ------------------------------------------------
        // determine the channel compression modes
        // ------------------------------------------------
        size_t dstCHA = CHA;
        if(downstream_coil_compression_num_modesKept_>0 && downstream_coil_compression_num_modesKept_<=CHA)
        {
            dstCHA = downstream_coil_compression_num_modesKept_;
        }
        else
        {
            std::vector<T> E(CHA, 0);
            long long cha;

            std::complex<T>* pRef = ref.begin();

            hoNDArray< std::complex<T> > dataCha;
            for (cha = 0; cha < (long long)CHA; cha++)
            {
                dataCha.create(RO, E1, pRef + cha*RO*E1);
                T v = Gadgetron::nrm2(dataCha);
                E[cha] = v*v;
            }

            for (cha = 1; cha < (long long)CHA; cha++)
            {
                if (std::abs(E[cha]) < downstream_coil_compression_thres_*std::abs(E[0]))
                {
                    break;
                }
            }

            dstCHA = cha;
        }

        GDEBUG_STREAM("--> perform_raw_data_recon, determine number of dst channels : " << dstCHA);

        ArrayType ref_src(RO, E1, CHA, ref.begin());
        ArrayType ref_dst(RO, E1, dstCHA, ref.begin());

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(ref_src, debug_folder_ + "raw_data_recon_ref_src" + suffix_);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(ref_dst, debug_folder_ + "raw_data_recon_ref_dst" + suffix_);

        // ------------------------------------------------
        // estimate coil map
        // ------------------------------------------------

        size_t start_RO = binning_obj_.sampling_.sampling_limits_[0].min_;
        size_t end_RO = binning_obj_.sampling_.sampling_limits_[0].max_;

        size_t start_E1 = binning_obj_.sampling_.sampling_limits_[1].min_;
        size_t end_E1 = binning_obj_.sampling_.sampling_limits_[1].max_;

        size_t lenRO = RO;
        if ((start_RO<RO) && (end_RO<RO) && (end_RO - start_RO + 1 < RO))
        {
            lenRO = (end_RO - start_RO + 1);
        }

        GDEBUG_STREAM("--> perform_raw_data_recon, RO sampling : [" << start_RO << " " << end_RO << "] out of " << RO);
        GDEBUG_STREAM("--> perform_raw_data_recon, RO sampling : [" << start_E1 << " " << end_E1 << "] out of " << E1);

        // compute ref coil map filter
        ArrayType filter_RO_ref_coi_map;
        Gadgetron::generate_symmetric_filter_ref(ref.get_size(0), start_RO, end_RO, filter_RO_ref_coi_map);

        ArrayType filter_E1_ref_coi_map;
        Gadgetron::generate_symmetric_filter_ref(ref.get_size(1), start_E1, end_E1, filter_E1_ref_coi_map);

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(filter_RO_ref_coi_map, debug_folder_ + "raw_data_recon_filter_RO_ref_coi_map" + suffix_);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(filter_E1_ref_coi_map, debug_folder_ + "raw_data_recon_filter_E1_ref_coi_map" + suffix_);

        // apply ref coil map filter
        ArrayType ref_coil_map_dst;
        Gadgetron::apply_kspace_filter_ROE1(ref_dst, filter_RO_ref_coi_map, filter_E1_ref_coi_map, ref_coil_map_dst);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(ref_coil_map_dst, debug_folder_ + "raw_data_recon_ref_coil_map_dst" + suffix_);

        // compute coil map
        ArrayType complex_im_coil_map;
        Gadgetron::hoNDFFT<T>::instance()->ifft2c(ref_coil_map_dst, complex_im_coil_map);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(complex_im_coil_map, debug_folder_ + "raw_data_recon_complex_im_coil_map" + suffix_);

        size_t ks = 7;
        size_t power = 3;

        if ( this->perform_timing_ ) { gt_timer_.start("--> perform_raw_data_recon, compute coil map"); }

        Gadgetron::coil_map_2d_Inati(complex_im_coil_map, coil_map, ks, power);

        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(coil_map, debug_folder_ + "raw_data_recon_Result_coil_map" + suffix_);

        // ------------------------------------------------
        // estimate grappa kernel
        // ------------------------------------------------

        size_t convKRO(1), convKE1(1), convKE2(1);

        std::vector<int> kE1, oE1;
        bool fitItself = true;
        Gadgetron::grappa2d_kerPattern(kE1, oE1, convKRO, convKE1, (size_t)binning_obj_.accel_factor_E1_, grappa_kSize_RO_, grappa_kSize_E1_, fitItself);

        ArrayType convKer(convKRO, convKE1, CHA, dstCHA);
        ArrayType kernelIm(RO, E1, CHA, dstCHA);

        Gadgetron::clear(convKer);
        Gadgetron::clear(kernelIm);

        if ( this->perform_timing_ ) { gt_timer_.start("--> perform_raw_data_recon, estimate grappa kernel"); }
        Gadgetron::grappa2d_calib_convolution_kernel(ref_src, ref_dst, (size_t)binning_obj_.accel_factor_E1_, grappa_reg_lamda_, grappa_kSize_RO_, grappa_kSize_E1_, convKer);
        Gadgetron::grappa2d_image_domain_kernel(convKer, RO, E1, kernelIm);
        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(convKer, debug_folder_ + "raw_data_recon_convKer" + suffix_);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kernelIm, debug_folder_ + "raw_data_recon_kernelIm" + suffix_);

        // compute unmixing coefficient
        ArrayType unmixing_coeff;
        hoNDArray<T> gFactor;
        Gadgetron::grappa2d_unmixing_coeff(kernelIm, coil_map, (size_t)binning_obj_.accel_factor_E1_, unmixing_coeff, gFactor);

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(unmixing_coeff, debug_folder_ + "raw_data_recon_unmixing_coeff" + suffix_);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array(gFactor, debug_folder_ + "raw_data_recon_gFactor" + suffix_);

        // ------------------------------------------------
        // perform reconstruction
        // ------------------------------------------------

        // compute aliased images
        ArrayType aliased_im(data);
        Gadgetron::hoNDFFT<T>::instance()->ifft2c(data, aliased_im);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(aliased_im, debug_folder_ + "raw_data_recon_aliased_im" + suffix_);

        // snr scaling
        float ROScalingFactor = (float)RO / (float)lenRO;
        float snr_scaling_ratio = (float)(std::sqrt(ROScalingFactor*binning_obj_.accel_factor_E1_));

        double grappaKernelCompensationFactor = 1.0 / (binning_obj_.accel_factor_E1_);
        Gadgetron::scal((T)(grappaKernelCompensationFactor*snr_scaling_ratio), aliased_im);

        GDEBUG_STREAM("--> perform_raw_data_recon, snr scaling factor : " << grappaKernelCompensationFactor*snr_scaling_ratio);

        // unwrapping
        if ( this->perform_timing_ ) { gt_timer_.start("--> perform_raw_data_recon, apply grappa kernel"); }

        Gadgetron::apply_unmix_coeff_aliased_image(aliased_im, unmixing_coeff, complex_image);

        if(this->use_multiple_channel_recon_)
        {
            // compute multi-channel full kspace
            Gadgetron::grappa2d_image_domain_unwrapping_aliased_image(aliased_im, kernelIm, full_kspace);
            Gadgetron::hoNDFFT<T>::instance()->fft2c(full_kspace);
        }
        else
        {
            // compute single-channel full kspace
            Gadgetron::hoNDFFT<T>::instance()->fft2c(complex_image, full_kspace);
        }

        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(complex_image, debug_folder_ + "raw_data_recon_Result_complex_image" + suffix_);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(full_kspace, debug_folder_ + "raw_data_recon_Result_full_kspace" + suffix_);
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::perform_raw_data_recon() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::estimate_time_stamps()
{
    try
    {
        hoNDArray< std::complex<T> >& kspace = binning_obj_.data_;
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t CHA = kspace.get_size(2);
        size_t N = kspace.get_size(3);
        size_t S = kspace.get_size(4);

        size_t startE1, endE1;
        startE1 = this->binning_obj_.sampling_.sampling_limits_[1].min_;
        endE1 = this->binning_obj_.sampling_.sampling_limits_[1].max_;

        size_t rE1 = endE1-startE1+1;

        hoNDArray< ISMRMRD::AcquisitionHeader >& header = binning_obj_.headers_;
        GADGET_CHECK_THROW(header.get_size(0)==E1);
        GADGET_CHECK_THROW(header.get_size(1)==N);
        GADGET_CHECK_THROW(header.get_size(2)==S);

        // initialize output arrays
        binning_obj_.ind_heart_beat_.create(E1, N, S);
        Gadgetron::clear(binning_obj_.ind_heart_beat_);

        binning_obj_.time_stamp_.create(E1, N, S);
        Gadgetron::fill(binning_obj_.time_stamp_, (float)-1);

        binning_obj_.cpt_time_stamp_.create(E1, N, S);
        Gadgetron::fill(binning_obj_.cpt_time_stamp_, (float)-1);

        // fill the time stamp
        size_t ii, e1, n, s;
        for (s=0; s<S; s++)
        {
            for (n=0; n<N; n++)
            {
                for (e1=0; e1<E1; e1++)
                {
                    if(header(e1,n,s).number_of_samples>0)
                    {
                        binning_obj_.time_stamp_(e1, n, s) = header(e1, n, s).acquisition_time_stamp * this->time_tick_;
                        binning_obj_.cpt_time_stamp_(e1, n, s) = header(e1, n, s).physiology_time_stamp[this->trigger_time_index_] * this->time_tick_;
                    }
                }
            }
        }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.time_stamp_, debug_folder_+"binning_obj_time_stamp" + suffix_);

        binning_obj_.cpt_time_ratio_.create(E1, N, S);
        Gadgetron::fill(binning_obj_.cpt_time_ratio_, (float)-1);

        binning_obj_.phs_time_stamp_.create(N, S);
        binning_obj_.phs_cpt_time_stamp_.create(N, S);
        binning_obj_.phs_cpt_time_ratio_.create(N, S);

        binning_obj_.starting_heart_beat_.resize(S);
        binning_obj_.ending_heart_beat_.resize(S);

        hoNDArray<long long> dummy(E1, N);

        for ( s=0; s<S; s++ )
        {
            std::stringstream os;
            os << "_S_" << s;

            // for current S
            hoNDArray<float> time_stamp(E1, N, binning_obj_.time_stamp_.begin()+s*E1*N);

            hoNDArray<float> cpt_time_stamp(E1, N, binning_obj_.cpt_time_stamp_.begin()+s*E1*N);
            hoNDArray<float> cpt_time_ratio(E1, N, binning_obj_.cpt_time_ratio_.begin()+s*E1*N);

            hoNDArray<int> ind_heart_beat(E1, N, binning_obj_.ind_heart_beat_.begin()+s*E1*N);

            HeartBeatIndexType startingHB, endingHB;

            hoNDArray<float> phs_time_stamp(N, 1, binning_obj_.phs_time_stamp_.begin()+s*N);
            hoNDArray<float> phs_cpt_time_stamp(N, 1, binning_obj_.phs_cpt_time_stamp_.begin()+s*N);
            hoNDArray<float> phs_cpt_time_ratio(N, 1, binning_obj_.phs_cpt_time_ratio_.begin()+s*N);

            // handle alternating sampling
            hoNDArray<float> cpt_time_stamp_flipped(E1, N);
            std::vector<bool> ascending;

            this->detect_and_flip_alternating_order(time_stamp, cpt_time_stamp, cpt_time_stamp_flipped, ascending);

            bool hasDescending = false;
            for ( n=0; n<N; n++ )
            {
                if ( ascending[n] == false )
                {
                    hasDescending = true;
                    break;
                }
            }

            if ( !debug_folder_.empty() ) gt_exporter_.export_array(time_stamp, debug_folder_+"time_stamp" + os.str() + suffix_);
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_stamp, debug_folder_+"cpt_time_stamp" + os.str() + suffix_);

            float meanRR(0);

            if ( !hasDescending )
            {
                if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_stamp, debug_folder_+"cpt_time_stamp" + os.str() + suffix_);

                this->process_time_stamps(time_stamp, cpt_time_stamp, 
                                cpt_time_ratio, phs_time_stamp, phs_cpt_time_stamp, phs_cpt_time_ratio, 
                                ind_heart_beat, startingHB, endingHB, meanRR);
            }
            else
            {
                GDEBUG_STREAM("CmrKSpaceBinning -- Alternating sampling strategy detected ... ");

                if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_stamp_flipped, debug_folder_+"cpt_time_stamp_flipped" + os.str() + suffix_);

                hoNDArray<float> cpt_time_ratio_tmp(E1, N);

                this->process_time_stamps(time_stamp, cpt_time_stamp_flipped, 
                                cpt_time_ratio_tmp, phs_time_stamp, phs_cpt_time_stamp, phs_cpt_time_ratio, 
                                ind_heart_beat, startingHB, endingHB, meanRR);

                if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_ratio_tmp, debug_folder_+"cpt_time_ratio_tmp" + os.str() + suffix_);

                cpt_time_stamp = cpt_time_stamp_flipped;
                cpt_time_ratio = cpt_time_ratio_tmp;

                for ( n=0; n<N; n++ )
                {
                    if ( !ascending[n] )
                    {
                        for ( e1=startE1; e1<=endE1; e1++ )
                        {
                            cpt_time_ratio(e1, n) = cpt_time_ratio_tmp(E1-1-e1, n);
                        }
                    }
                }
            }

            if ( !debug_folder_.empty() ) gt_exporter_.export_array(time_stamp, debug_folder_+"time_stamp_after_processing" + os.str() + suffix_);
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_stamp, debug_folder_+"cpt_time_stamp_after_processing" + os.str() + suffix_);
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_ratio, debug_folder_+"cpt_time_ratio" + os.str() + suffix_);
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(phs_time_stamp, debug_folder_+"phs_time_stamp" + os.str() + suffix_);
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(phs_cpt_time_stamp, debug_folder_+"phs_cpt_time_stamp" + os.str() + suffix_);
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(phs_cpt_time_ratio, debug_folder_+"phs_cpt_time_ratio" + os.str() + suffix_);

            binning_obj_.starting_heart_beat_[s] = startingHB;
            binning_obj_.ending_heart_beat_[s] = endingHB;
        }

        // compute desired cardiac phase time stamp and time ratio
        size_t total_num_HB = 0;
        std::vector<float> starting_HB_time_stamp;
        std::vector<float> ending_HB_time_stamp;

        size_t ind;
        for ( s=0; s<S; s++ )
        {
            size_t num_HB_S = binning_obj_.starting_heart_beat_[s].size();

            for ( ind=1; ind<num_HB_S-1; ind++ )
            {
                size_t start_e1 = binning_obj_.starting_heart_beat_[s][ind].first;
                size_t end_e1 = binning_obj_.ending_heart_beat_[s][ind].first;

                size_t start_n = binning_obj_.starting_heart_beat_[s][ind].second;
                size_t end_n = binning_obj_.ending_heart_beat_[s][ind].second;

                size_t start_ind = start_e1 + start_n * E1;
                size_t end_ind = end_e1 + end_n * E1;

                float start_cpt = -1;
                for (ii=start_ind; ii<=end_ind; ii++)
                {
                    if(binning_obj_.cpt_time_stamp_(ii)>=0)
                    {
                        start_cpt = binning_obj_.cpt_time_stamp_(ii);
                        break;
                    }
                }

                float end_cpt = -1;
                for (ii=end_ind; ii>=start_ind; ii--)
                {
                    if(binning_obj_.cpt_time_stamp_(ii)>=0)
                    {
                        end_cpt = binning_obj_.cpt_time_stamp_(ii);
                        break;
                    }
                }

                starting_HB_time_stamp.push_back(start_cpt);
                ending_HB_time_stamp.push_back(end_cpt);
                total_num_HB++;
            }
        }

        if ( total_num_HB == 0 ) total_num_HB = 1;

        // compute mean RR
        std::sort(starting_HB_time_stamp.begin(), starting_HB_time_stamp.end());
        float mean_starting_HB_time_stamp = starting_HB_time_stamp[starting_HB_time_stamp.size()/2];

        std::sort(ending_HB_time_stamp.begin(), ending_HB_time_stamp.end());
        float mean_ending_HB_time_stamp = ending_HB_time_stamp[ending_HB_time_stamp.size()/2];

        binning_obj_.mean_RR_ = mean_ending_HB_time_stamp - mean_starting_HB_time_stamp;
        GDEBUG_STREAM("CmrKSpaceBinning -- mean RR : " << binning_obj_.mean_RR_ << "ms -- start and end heart beat time stamp : [" << mean_starting_HB_time_stamp << " " << mean_ending_HB_time_stamp << "]");

        // compute output cardiac phase time stamp and time ratio
        float step_size_time_stamp = binning_obj_.mean_RR_ / binning_obj_.output_N_;
        float step_size_time_ratio = 1.0f / binning_obj_.output_N_;

        binning_obj_.desired_cardiac_phases_.resize(binning_obj_.output_N_, 0);
        binning_obj_.desired_cpt_.resize(binning_obj_.output_N_, 0);

        for ( ind=0; ind<binning_obj_.output_N_; ind++ )
        {
            if( ind == 0)
            {
                binning_obj_.desired_cardiac_phases_[0] = step_size_time_stamp/2;
                binning_obj_.desired_cpt_[0] = step_size_time_ratio/2;
            }
            else
            {
                binning_obj_.desired_cardiac_phases_[ind] = binning_obj_.desired_cardiac_phases_[ind-1] + step_size_time_stamp;
                binning_obj_.desired_cpt_[ind] = binning_obj_.desired_cpt_[ind-1] + step_size_time_ratio;
            }

            GDEBUG_STREAM("CmrKSpaceBinning -- desired_cardiac_phases -- " << ind << " - " << binning_obj_.desired_cardiac_phases_[ind]);
            GDEBUG_STREAM("CmrKSpaceBinning -- desired_cpt            -- " << ind << " - " << binning_obj_.desired_cpt_[ind]);
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::estimate_time_stamps() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::estimate_respiratory_navigator()
{
    try
    {
        // estimate respiratory signal from binning_obj_.complex_image_raw_
        ArrayType& complex_image = binning_obj_.complex_image_raw_;

        size_t RO = complex_image.get_size(0);
        size_t E1 = complex_image.get_size(1);
        size_t CHA = complex_image.get_size(2);
        size_t N = complex_image.get_size(3);
        size_t S = complex_image.get_size(4);

        binning_obj_.navigator_.create(E1, N, S);
        Gadgetron::clear(binning_obj_.navigator_);

        hoNDArray<float>& navigator = binning_obj_.navigator_;

        NavigatorRoiType roi;
        roi.resize(S);

        Gadgetron::hoImageRegContainer2DRegistration<ImageType, ImageType, double> reg;

        size_t s, n, e1;
        for (s=0; s<S; s++)
        {
            std::stringstream os;
            os << "_S_" << s;

            ArrayType im(RO, E1, N, complex_image.begin()+s*RO*E1*N);
            hoNDArray<float> navigator_s(E1, N, navigator.begin()+s*E1*N);

            // get the magnitude of image
            hoNDArray<T> mag;
            Gadgetron::abs(im, mag);

            if ( !debug_folder_.empty() ) gt_exporter_.export_array(mag, debug_folder_+"navigator_mag" + os.str() + suffix_);

            // detect key frame
            size_t key_frame;
            Gadgetron::find_key_frame_SSD_2DT(mag, key_frame);
            GDEBUG_STREAM("Find key frame " << key_frame << " for S " << s);

            // perform motion correction
            bool bidirectional_moco = false;
            bool warp_images = false;
            if ( !debug_folder_.empty() ) warp_images = true;

            Gadgetron:: perform_moco_fixed_key_frame_2DT(mag, key_frame, respiratory_navigator_moco_reg_strength_, respiratory_navigator_moco_iters_, bidirectional_moco, warp_images, reg);

            // get the moco results
            hoNDArray<double> dx, dy;
            reg.deformation_field_[0].to_NDArray(0, dx);
            reg.deformation_field_[1].to_NDArray(0, dy);

            if ( !debug_folder_.empty() ) gt_exporter_.export_array(dx, debug_folder_+"navigator_mag_moco_dx" + os.str() + suffix_);
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(dy, debug_folder_+"navigator_mag_moco_dy" + os.str() + suffix_);

            if ( !debug_folder_.empty() )
            {
                hoNDArray<T> magMoCo(mag);
                reg.warped_container_.to_NDArray(0, magMoCo);
                gt_exporter_.export_array(magMoCo, debug_folder_+"navigator_mag_moco" + os.str() + suffix_);
            }

            // select the best region and get the navigator signal

            // generate the region
            std::vector<long long> starting_pos_RO, starting_pos_E1;

            long long patch_size_RO = respiratory_navigator_patch_size_RO_;
            long long patch_size_E1 = respiratory_navigator_patch_size_E1_;

            long long patch_step_size_RO = respiratory_navigator_patch_step_size_RO_;
            long long patch_step_size_E1 = respiratory_navigator_patch_step_size_E1_;

            long long pos = 0;
            while ( pos <= RO-patch_size_RO )
            {
                starting_pos_RO.push_back(pos);
                pos += patch_step_size_RO;
            }

            pos = 0;
            while ( pos <= E1-patch_size_E1 )
            {
                starting_pos_E1.push_back(pos);
                pos += patch_step_size_E1;
            }

            unsigned long long num_RO = starting_pos_RO.size();
            unsigned long long num_E1 = starting_pos_E1.size();

            // compute the mean deformation within the region
            hoNDArray<float> mean_deform_RO(num_RO, num_E1, N);
            hoNDArray<float> mean_deform_E1(num_RO, num_E1, N);
            hoNDArray<double> deform_patch(patch_size_RO, patch_size_E1);

            size_t l, c;
            for ( n=0; n<N; n++ )
            {
                for ( l=0; l<num_E1; l++ )
                {
                    for ( c=0; c<num_RO; c++ )
                    {
                        long long pl, k;

                        // X
                        // copy deformation to deform_patch
                        for ( pl=0; pl<patch_size_E1; pl++ )
                        {
                            long long offset = dx.calculate_offset( starting_pos_RO[c], starting_pos_E1[l]+pl, n );
                            memcpy(deform_patch.begin()+pl*patch_size_RO, dx.begin()+offset, sizeof(double)*patch_size_RO);
                        }

                        // compute the mean deformation
                        float totalDeform = 0.0f;
                        for ( k=0; k<patch_size_RO*patch_size_E1; k++ )
                        {
                            totalDeform += deform_patch(k);
                        }

                        mean_deform_RO(c, l, n) = totalDeform/(patch_size_RO*patch_size_E1);

                        // Y
                        // copy deformation to deform_patch
                        for ( pl=0; pl<patch_size_E1; pl++ )
                        {
                            long long offset = dy.calculate_offset( starting_pos_RO[c], starting_pos_E1[l]+pl, n );
                            memcpy(deform_patch.begin()+pl*patch_size_RO, dy.begin()+offset, sizeof(double)*patch_size_RO);
                        }

                        // compute the mean deformation
                        totalDeform = 0.0f;
                        for ( k=0; k<patch_size_RO*patch_size_E1; k++ )
                        {
                            totalDeform += deform_patch(k);
                        }

                        mean_deform_E1(c, l, n) = totalDeform/(patch_size_RO*patch_size_E1);
                    }
                }
            }

            // at this moment, do it simply by computing the variance for every region and pick the one with largest variance
            ho2DArray<float> var_mean_deform_RO(num_RO, num_E1);
            ho2DArray<float> var_mean_deform_E1(num_RO, num_E1);
            for ( l=0; l<num_E1; l++ )
            {
                for ( c=0; c<num_RO; c++ )
                {
                    hoNDArray<float> buf(N);

                    for ( n=0; n<N; n++ )
                    {
                        buf(n) = mean_deform_RO(c, l, n);
                    }

                    var_mean_deform_RO(c, l) = Gadgetron::stddev(&buf);

                    for ( n=0; n<N; n++ )
                    {
                        buf(n) = mean_deform_E1(c, l, n);
                    }
                    var_mean_deform_E1(c, l) = Gadgetron::stddev(&buf);
                }
            }

            // find the region with maximal variance
            size_t max_var_ind_RO, max_var_ind_E1;
            float max_var_RO, max_var_E1;
            Gadgetron::maxAbsolute(var_mean_deform_RO, max_var_RO, max_var_ind_RO);
            Gadgetron::maxAbsolute(var_mean_deform_E1, max_var_E1, max_var_ind_E1);

            std::vector<size_t> regionInd = var_mean_deform_RO.calculate_index(max_var_ind_RO);
            if ( max_var_E1 > max_var_RO )
            {
                regionInd = var_mean_deform_E1.calculate_index(max_var_ind_E1);
            }

            roi[s].first[0] = starting_pos_RO[regionInd[0]];
            roi[s].first[1] = starting_pos_E1[regionInd[1]];

            roi[s].second[0] = starting_pos_RO[regionInd[0]]+patch_size_RO;
            roi[s].second[1] = starting_pos_E1[regionInd[1]]+patch_size_E1;

            GDEBUG_CONDITION_STREAM(this->verbose_, "Respiratory ROI - RO - E1 [" 
                << roi[s].first(0) << "," << roi[s].first(1) << "]; [" << roi[s].second(0) << " -- " << roi[s].second(1) << "]");

            for ( n=0; n<N; n++ )
            {
                float deformValue = mean_deform_RO(regionInd[0], regionInd[1], n);
                if ( max_var_E1 > max_var_RO )
                {
                    deformValue = mean_deform_E1(regionInd[0], regionInd[1], n);
                }

                for ( e1=0; e1<E1; e1++ )
                {
                    navigator_s(e1, n) =  deformValue;
                }
                GDEBUG_CONDITION_STREAM(verbose_, "Respirator navigator - n " << n << ", s " << s << " : " << deformValue);
            }

            if ( !debug_folder_.empty() ) gt_exporter_.export_array(navigator_s, debug_folder_ + "navigator"  + os.str() + suffix_);
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::estimate_respiratory_navigator() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::find_best_heart_beat(std::vector<size_t>& bestHB)
{
    try
    {
        ArrayType& complex_image = binning_obj_.complex_image_raw_;

        size_t RO = complex_image.get_size(0);
        size_t E1 = complex_image.get_size(1);
        size_t CHA = complex_image.get_size(2);
        size_t N = complex_image.get_size(3);
        size_t S = complex_image.get_size(4);

        hoNDArray<float>& navigator = binning_obj_.navigator_;
        hoNDArray<float>& cpt_time_stamp = binning_obj_.cpt_time_stamp_;

        bestHB.resize(S);

        size_t s, n, e1;
        for (s=0; s<S; s++)
        {
            std::stringstream os;
            os << "_S_" << s;

            hoNDArray<float> navigator_s(E1, N, navigator.begin()+s*E1*N);
            hoNDArray<float> cpt_time_stamp_s(E1, N, cpt_time_stamp.begin()+s*E1*N);

            size_t numOfHB = binning_obj_.starting_heart_beat_[s].size();
            bestHB[s] = 0;
            float min_var_resp = (float)(1e15);

            size_t HB;
            for ( HB=1; HB<numOfHB-1; HB++ ) // do not pick the first the last heart beat
            {
                // reject irregular heart beat
                float RRInterval;
                this->compute_RRInterval(s, HB, RRInterval);

                if ( RRInterval < binning_obj_.mean_RR_ * (1.0-this->arrhythmia_rejector_factor_) )
                    continue;

                if ( RRInterval > binning_obj_.mean_RR_ * (1.0+this->arrhythmia_rejector_factor_) )
                    continue;

                float mean_resp, var_resp, min_resp, max_resp;
                this->compute_metrics_navigator_heart_beat(s, HB, mean_resp, var_resp, min_resp, max_resp);

                if ( var_resp < min_var_resp )
                {
                    min_var_resp = var_resp;
                    bestHB[s] = HB;
                }
            }

            GDEBUG_CONDITION_STREAM(this->verbose_, "find best heart beat : " << bestHB[s] << " for S " << s);
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::find_best_heart_beat() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::reject_irregular_heart_beat()
{
    try
    {
        ArrayType& complex_image = binning_obj_.complex_image_raw_;

        size_t RO = complex_image.get_size(0);
        size_t E1 = complex_image.get_size(1);
        size_t CHA = complex_image.get_size(2);
        size_t N = complex_image.get_size(3);
        size_t S = complex_image.get_size(4);

        hoNDArray<float>& cpt_time_stamp = binning_obj_.cpt_time_stamp_;
        hoNDArray<float>& cpt_time_ratio = binning_obj_.cpt_time_ratio_;

        size_t s, n, e1;
        for (s=0; s<S; s++)
        {
            std::stringstream os;
            os << "_S_" << s;

            hoNDArray<float> cpt_time_stamp_s(E1, N, cpt_time_stamp.begin()+s*E1*N);
            hoNDArray<float> cpt_time_ratio_s(E1, N, cpt_time_ratio.begin()+s*E1*N);

            size_t numOfHB = binning_obj_.starting_heart_beat_[s].size();

            size_t HB;
            for ( HB=0; HB<numOfHB-1; HB++ ) // do not pick the first the last heart beat
            {
                float RRInterval;
                this->compute_RRInterval(s, HB, RRInterval);
                GDEBUG_CONDITION_STREAM(this->verbose_, "Heart beat " << HB << " - S " << s << " - RRInterval - " << RRInterval);

                bool reject = false;
                if ( RRInterval < binning_obj_.mean_RR_ * (1.0-this->arrhythmia_rejector_factor_) )
                    reject = true;

                if ( RRInterval > binning_obj_.mean_RR_ * (1.0+this->arrhythmia_rejector_factor_) )
                    reject = true;

                if(reject)
                {
                    GDEBUG_CONDITION_STREAM(this->verbose_, "Heart beat " << HB << " is rejected due to RR interval - " << RRInterval << " - meanRR " << binning_obj_.mean_RR_);

                    size_t sInd = binning_obj_.starting_heart_beat_[s][HB].second*E1 + binning_obj_.starting_heart_beat_[s][HB].first;
                    size_t eInd = binning_obj_.ending_heart_beat_[s][HB].second*E1 + binning_obj_.ending_heart_beat_[s][HB].first;

                    size_t ind;
                    for ( ind=sInd; ind<=eInd; ind++ )
                    {
                        cpt_time_stamp_s(ind) = (float)(-1e3);
                        cpt_time_ratio_s(ind) = (float)(-1e3);
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::reject_irregular_heart_beat() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::detect_and_flip_alternating_order(const hoNDArray<float>& time_stamp, const hoNDArray<float>& cpt_time_stamp, hoNDArray<float>& cpt_time_stamp_flipped, std::vector<bool>& ascending)
{
    try
    {
        long long E1 = cpt_time_stamp.get_size(0);
        long long N = cpt_time_stamp.get_size(1);

        ascending.resize(N, true);

        long long n, e1;

        for ( n=0; n<N; n++ )
        {
            float v_first = -1;
            for ( e1=0; e1<E1; e1++ )
            {
                if ( time_stamp(e1, n) > -1 )
                {
                    v_first = time_stamp(e1, n);
                    break;
                }
            }

            float v_last = -1;
            for ( e1=E1-1; e1>=0; e1-- )
            {
                if ( time_stamp(e1, n) > -1 )
                {
                    v_last = time_stamp(e1, n);
                    break;
                }
            }

            if ( v_last < v_first ) ascending[n] = false;
        }

        cpt_time_stamp_flipped = cpt_time_stamp;
        Gadgetron::fill(cpt_time_stamp_flipped, (float)(-1) );

        for ( n=0; n<N; n++ )
        {
            if ( !ascending[n] )
            {
                for ( e1=E1-1; e1>=0; e1-- )
                {
                    if ( cpt_time_stamp(e1, n) > -1 )
                    {
                        cpt_time_stamp_flipped(E1-1-e1, n) = cpt_time_stamp(e1, n);
                    }
                }
            }
            else
            {
                for ( e1=0; e1<E1; e1++ )
                {
                    if ( cpt_time_stamp(e1, n) > -1 )
                    {
                        cpt_time_stamp_flipped(e1, n) = cpt_time_stamp(e1, n);
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::detect_and_flip_alternating_order() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::process_time_stamps(hoNDArray<float>& time_stamp, hoNDArray<float>& cpt_time_stamp, 
                                            hoNDArray<float>& cpt_time_ratio, hoNDArray<float>& phs_time_stamp, 
                                            hoNDArray<float>& phs_cpt_time_stamp, hoNDArray<float>& phs_cpt_time_ratio, 
                                            hoNDArray<int>& indHeartBeat, HeartBeatIndexType& startingHB, 
                                            HeartBeatIndexType& endingHB, float& meanRRInterval)
{
    try
    {
        size_t E1 = time_stamp.get_size(0);
        size_t N = time_stamp.get_size(1);

        size_t startE1, endE1;
        startE1 = this->binning_obj_.sampling_.sampling_limits_[1].min_;
        endE1 = this->binning_obj_.sampling_.sampling_limits_[1].max_;

        size_t rE1 = endE1-startE1+1;

        size_t e1, n, ind, ii;

        startingHB.resize(1);
        endingHB.resize(1);

        // --------------------------------------------------------
        // fit a line along time stamp
        // --------------------------------------------------------

        size_t num_acq_read_outs = 0;
        for ( n=0; n<N; n++ )
        {
            for ( e1=0; e1<E1; e1++ )
            {
                if ( time_stamp(e1, n) > 0 )
                {
                    num_acq_read_outs++;
                }
            }
        }

        GDEBUG_CONDITION_STREAM(this->verbose_, " Number of acquired lines : " << num_acq_read_outs);

        Gadgetron::correct_time_stamp_with_fitting(time_stamp, startE1, endE1);

        if ( !debug_folder_.empty() ) gt_exporter_.export_array(time_stamp, debug_folder_+"time_stamp_after_fitting" + suffix_);

        // --------------------------------------------------------
        // cpt time stamps
        // --------------------------------------------------------
        std::vector<size_t> start_e1_hb, end_e1_hb, start_n_hb, end_n_hb;
        Gadgetron::detect_heart_beat_with_time_stamp(cpt_time_stamp, indHeartBeat, start_e1_hb, end_e1_hb, start_n_hb, end_n_hb);

        size_t numOfHB = start_e1_hb.size();
        GADGET_CHECK_THROW( numOfHB > 1 );

        // --------------------------------------------------------
        // fill the starting and endingHB
        // --------------------------------------------------------
        startingHB.resize(numOfHB);
        endingHB.resize(numOfHB);

        std::vector<size_t> start, end;
        for ( ind=0; ind<numOfHB; ind++ )
        {
            startingHB[ind].first = start_e1_hb[ind];
            startingHB[ind].second = start_n_hb[ind];

            endingHB[ind].first = end_e1_hb[ind];
            endingHB[ind].second = end_n_hb[ind];

            GDEBUG_CONDITION_STREAM(verbose_, "Heart beat " << ind 
                << " - start " << " [" << startingHB[ind].first << " " << startingHB[ind].second 
                << "] - end  " << " [" << endingHB[ind].first << " " << endingHB[ind].second << "]");
        }

        // --------------------------------------------------------
        // correct cpt time stamp
        // --------------------------------------------------------
        Gadgetron::correct_heart_beat_time_stamp_with_fitting(cpt_time_stamp, indHeartBeat, startE1, endE1, start_e1_hb, end_e1_hb, start_n_hb, end_n_hb);

        // --------------------------------------------------------
        // fill per phase time stamp  and phase cpt time stamp
        // --------------------------------------------------------
        Gadgetron::compute_phase_time_stamp(time_stamp, cpt_time_stamp, startE1, endE1, phs_time_stamp, phs_cpt_time_stamp);

        // --------------------------------------------------------
        // compute the mean maximal cpt
        // --------------------------------------------------------
        std::vector<size_t> ind_HB_start(numOfHB);
        std::vector<size_t> ind_HB_end(numOfHB);

        for ( ind=0; ind<numOfHB; ind++ )
        {
            ind_HB_start[ind] = start_e1_hb[ind] + start_n_hb[ind] * E1;
            ind_HB_end[ind] = end_e1_hb[ind] + end_n_hb[ind] * E1;
        }

        float mean_max_cpt = 0;
        for ( ind=0; ind<numOfHB-1; ind++ )
        {
            // find the max and min cpt time for this heart beat
            float maxCPT = -1e10f;
            for ( ii=ind_HB_start[ind]; ii<=ind_HB_end[ind]; ii++ )
            {
                size_t n = ii / E1;
                size_t e1 = ii - n*E1;

                float v = cpt_time_stamp(e1, n);
                if ( v > -1 )
                {
                    if ( v > maxCPT )
                        maxCPT = v;
                }
            }

            mean_max_cpt += maxCPT;
        }
        mean_max_cpt /= (numOfHB-1);

        // --------------------------------------------------------
        // fill the cpt time ratio
        // --------------------------------------------------------
        for ( ind=0; ind<numOfHB; ind++ )
        {
            // find the max and min cpt time for this heart beat
            float maxCPT = -1e10f;
            float minCPT = 1e10f;
            for ( ii=ind_HB_start[ind]; ii<=ind_HB_end[ind]; ii++ )
            {
                size_t n = ii / E1;
                size_t e1 = ii - n*E1;

                float v = cpt_time_stamp(e1, n);

                if ( v > -1 )
                {
                    if ( v > maxCPT )
                        maxCPT = v;

                    if ( v < minCPT )
                        minCPT = v;
                }
            }

            if ( std::abs(maxCPT-minCPT) < FLT_EPSILON )
            {
                maxCPT = minCPT + 1.0f;
            }

            // handle the last heart beat
            if ( ind == numOfHB-1 )
            {
                if ( maxCPT < mean_max_cpt )
                {
                    maxCPT = mean_max_cpt;
                }
            }

            // fill in the cptTimeRatio
            for ( ii=ind_HB_start[ind]; ii<=ind_HB_end[ind]; ii++ )
            {
                size_t n = ii / E1;
                size_t e1 = ii - n*E1;

                if(e1>=startE1 && e1<=endE1)
                {
                    cpt_time_ratio(e1, n) = (cpt_time_stamp(e1, n) - minCPT) / (maxCPT-minCPT);
                }
            }
        }

        // --------------------------------------------------------
        // fill the phase cpt time ratio
        // --------------------------------------------------------
        for ( n=0; n<N; n++ )
        {
            std::vector<float> buf(rE1);
            for ( e1=startE1; e1<=endE1; e1++ )
            {
                buf[e1-startE1] = cpt_time_ratio(e1, n);
            }

            std::sort(buf.begin(), buf.end());
            phs_cpt_time_ratio(n, 0) = (float)(0.5*(buf[E1/2-startE1-1]+buf[E1/2-startE1]));
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::process_time_stamps() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::compute_RRInterval(size_t s, size_t HB, float& RRInterval )
{
    try
    {
        GADGET_CHECK_THROW(HB<binning_obj_.starting_heart_beat_[s].size());

        size_t E1 = binning_obj_.cpt_time_stamp_.get_size(0);
        size_t N = binning_obj_.cpt_time_stamp_.get_size(1);

        size_t startOffset = binning_obj_.starting_heart_beat_[s][HB].second*E1 + binning_obj_.starting_heart_beat_[s][HB].first;
        size_t endOffset = binning_obj_.ending_heart_beat_[s][HB].second*E1 + binning_obj_.ending_heart_beat_[s][HB].first;

        float min_time_stamp = 1e12f;
        float max_time_stamp = -1;

        size_t ind;
        for ( ind=startOffset; ind<=endOffset; ind++ )
        {
            float v = binning_obj_.cpt_time_stamp_.at(ind);
            if ( v>=0 && v < min_time_stamp ) min_time_stamp = v;
            if ( v>=0 && v > max_time_stamp ) max_time_stamp = v;
        }

        RRInterval = max_time_stamp - min_time_stamp;
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::compute_RRInterval() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::compute_metrics_navigator_heart_beat(size_t s, size_t HB, float& mean_resp, float& var_resp, float& min_resp, float& max_resp)
{
    try
    {
        GADGET_CHECK_THROW(HB<binning_obj_.starting_heart_beat_[s].size());

        size_t E1 = binning_obj_.navigator_.get_size(0);
        size_t N = binning_obj_.navigator_.get_size(1);

        size_t sInd = binning_obj_.starting_heart_beat_[s][HB].second*E1 + binning_obj_.starting_heart_beat_[s][HB].first;
        size_t eInd = binning_obj_.ending_heart_beat_[s][HB].second*E1 + binning_obj_.ending_heart_beat_[s][HB].first;

        hoNDArray<float> nav(eInd-sInd+1);

        size_t ind;
        for ( ind=sInd; ind<=eInd; ind++ )
        {
            nav(ind-sInd) = binning_obj_.navigator_(ind);
        }

        var_resp = Gadgetron::stddev(&nav);
        mean_resp = Gadgetron::mean(&nav);
        min_resp = Gadgetron::min(&nav);
        max_resp = Gadgetron::max(&nav);
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::compute_metrics_navigator_heart_beat() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::compute_kspace_binning(const std::vector<size_t>& bestHB, std::vector<size_t>& slices_not_processing)
{
    try
    {
        ArrayType& full_kspace_raw = binning_obj_.full_kspace_raw_;
        ArrayType& coil_map_raw = binning_obj_.coil_map_raw_;
        ArrayType& complex_image_raw = binning_obj_.complex_image_raw_;

        size_t RO = full_kspace_raw.get_size(0);
        size_t E1 = full_kspace_raw.get_size(1);
        size_t CHA = full_kspace_raw.get_size(2);
        size_t N = full_kspace_raw.get_size(3);
        size_t S = full_kspace_raw.get_size(4);

        hoNDArray<T> mag;
        Gadgetron::abs(complex_image_raw, mag);

        if ( !debug_folder_.empty() ) gt_exporter_.export_array(mag, debug_folder_ + "complex_image_raw_mag" + suffix_);

        size_t dstN = binning_obj_.output_N_;
        binning_obj_.kspace_binning_.create(RO, E1, CHA, binning_obj_.output_N_, S);
        Gadgetron::clear(binning_obj_.kspace_binning_);

        binning_obj_.kspace_binning_wider_.create(RO, E1, CHA, binning_obj_.output_N_, S);
        Gadgetron::clear(binning_obj_.kspace_binning_wider_);

        binning_obj_.kspace_binning_image_domain_average_.create(RO, E1, CHA, binning_obj_.output_N_, S);
        Gadgetron::clear(binning_obj_.kspace_binning_image_domain_average_);

        binning_obj_.kspace_binning_hit_count_.create(E1, dstN, S);
        Gadgetron::clear(binning_obj_.kspace_binning_hit_count_);

        slices_not_processing.clear();

        size_t e1, n, s, ii;

        if(this->kspace_binning_interpolate_heart_beat_images_)
        {
            for (s=0; s<S; s++)
            {
                std::stringstream os;
                os << "_S_" << s;

                size_t bestHB_S = bestHB[s];

                size_t startE1 = binning_obj_.starting_heart_beat_[s][bestHB_S].first;
                size_t startN = binning_obj_.starting_heart_beat_[s][bestHB_S].second;

                size_t endE1 = binning_obj_.ending_heart_beat_[s][bestHB_S].first;
                size_t endN = binning_obj_.ending_heart_beat_[s][bestHB_S].second;

                size_t ori_endN = endN;
                size_t ori_startN = endN;

                // avoid cross the RR wav
                if ( binning_obj_.phs_cpt_time_stamp_(startN, s) > binning_obj_.phs_cpt_time_stamp_(startN+1, s) )
                {
                    startN++;
                }

                if ( binning_obj_.phs_cpt_time_stamp_(endN, s) < binning_obj_.phs_cpt_time_stamp_(endN-1, s) )
                {
                    endN--;
                }

                if ( endN <= startN )
                {
                    GERROR_STREAM("KSpace binning for S " << s << " - endN <= startN - " << endN << " <= " << startN );
                    GWARN_STREAM("Please consider to reduce temporal footprint of raw image series, if heart rate is too high ... " );

                    endN = ori_endN;
                    startN = ori_startN;

                    if ( endN <= startN )
                    {
                        GERROR_STREAM("KSpace binning for S " << s << " - endN <= startN - " << endN << " <= " << startN );
                        GERROR_STREAM("Slice " << s << " will not be processed ... ");
                        slices_not_processing.push_back(s);
                        continue;
                    }
                }

                // ----------------------------------------
                // get the best HB mag images
                // ----------------------------------------
                size_t num_images_bestHB = endN - startN + 1;
                hoNDArray<T> mag_bestHB(RO, E1, num_images_bestHB, mag.begin()+s*RO*E1*N+startN*RO*E1);

                if ( !debug_folder_.empty() ) gt_exporter_.export_array(mag_bestHB, debug_folder_ + "mag_bestHB" + os.str() + suffix_);

                // ----------------------------------------
                // get the respiratory location of the best HB
                // ----------------------------------------
                float mean_resp_bestHB(0), var_resp_bestHB(0), min_resp_bestHB(0), max_resp_bestHB(0);
                this->compute_metrics_navigator_heart_beat(s, bestHB[s], mean_resp_bestHB, var_resp_bestHB, min_resp_bestHB, max_resp_bestHB);

                float min_resp(0), max_resp(0);

                hoNDArray<float> navigator_S(E1, N, binning_obj_.navigator_.begin()+s*E1*N);
                min_resp = Gadgetron::min(&navigator_S);
                max_resp = Gadgetron::max(&navigator_S);

                // ----------------------------------------
                // interpolate best HB images to desired cardiac time ratio
                // ----------------------------------------
                std::vector<float> cpt_time_ratio_bestHB(num_images_bestHB);
                for ( n=startN; n<=endN; n++ )
                {
                    cpt_time_ratio_bestHB[n-startN] = binning_obj_.phs_cpt_time_ratio_(n, s);
                }

                hoNDArray<T> mag_bestHB_at_desired_cpt;
                this->interpolate_best_HB_images(cpt_time_ratio_bestHB, mag_bestHB, binning_obj_.desired_cpt_, mag_bestHB_at_desired_cpt);

                if ( !debug_folder_.empty() ) gt_exporter_.export_array(mag_bestHB_at_desired_cpt, debug_folder_ + "mag_bestHB_at_desired_cpt" + os.str() + suffix_);

                // ----------------------------------------
                // compute acceptance range
                // ----------------------------------------
                float accepted_nav_wider[2];

                float wider_nav_window = 1.3*this->kspace_binning_navigator_acceptance_window_;
                if (wider_nav_window>= 0.85) wider_nav_window = 0.85;

                accepted_nav_wider[0] = (float)(mean_resp_bestHB - 0.5*wider_nav_window*(max_resp-min_resp));
                accepted_nav_wider[1] = (float)(mean_resp_bestHB + 0.5*wider_nav_window*(max_resp-min_resp));

                if ( accepted_nav_wider[0] < min_resp )
                {
                    float delta = min_resp-accepted_nav_wider[0];
                    accepted_nav_wider[0] += delta;
                    accepted_nav_wider[1] += delta;
                }

                if ( accepted_nav_wider[1] > max_resp )
                {
                    float delta = accepted_nav_wider[1] - max_resp;
                    accepted_nav_wider[0] -= delta;
                    accepted_nav_wider[1] -= delta;
                }

                // ------------------------------------------------

                float accepted_nav[2];
                accepted_nav[0] = (float)(mean_resp_bestHB - 0.5*this->kspace_binning_navigator_acceptance_window_*(max_resp-min_resp));
                accepted_nav[1] = (float)(mean_resp_bestHB + 0.5*this->kspace_binning_navigator_acceptance_window_*(max_resp-min_resp));

                if ( accepted_nav[0] < min_resp )
                {
                    float delta = min_resp-accepted_nav[0];
                    accepted_nav[0] += delta;
                    accepted_nav[1] += delta;
                }

                if ( accepted_nav[1] > max_resp )
                {
                    float delta = accepted_nav[1] - max_resp;
                    accepted_nav[0] -= delta;
                    accepted_nav[1] -= delta;
                }

                // ----------------------------------------
                // for every destination N, compute images falling into its bin
                // ----------------------------------------
                std::vector < std::vector<size_t> > selected_images_wider(dstN);
                std::vector < std::vector<size_t> > selected_images(dstN);
                for (n=0; n<dstN; n++)
                {
                    float desired_cpt_time_ratio = binning_obj_.desired_cpt_[n];

                    std::vector<size_t> selected;

                    this->select_images_with_navigator(desired_cpt_time_ratio, accepted_nav_wider, n, s, selected);
                    selected_images_wider[n] = selected;

                    // --------------------------

                    size_t num = selected_images_wider[n].size();
                    for (size_t jj=0; jj<num; jj++)
                    {
                        float nav = binning_obj_.navigator_(E1/2, selected_images_wider[n][jj], s);

                        if(nav>=accepted_nav[0] && nav<=accepted_nav[1])
                        {
                            selected_images[n].push_back(selected_images_wider[n][jj]);
                        }
                    }

                    GDEBUG_CONDITION_STREAM(this->verbose_, "num of images selected for [n, s] : [" << n << ", " << s << "] is " << selected_images[n].size() << " out of " << selected_images_wider[n].size());
                }

                // ----------------------------------------
                // perform motion correction to compute deformation fields
                // ----------------------------------------
                hoNDArray<T> mag_s(RO, E1, N, mag.begin()+s*RO*E1*N);
                DeformationFieldContinerType deform[2];

                if ( this->perform_timing_ ) { gt_timer_local_.start("perform_moco_selected_images_with_best_heart_beat ... "); }
                this->perform_moco_selected_images_with_best_heart_beat(selected_images_wider, mag_s, mag_bestHB_at_desired_cpt, deform[0], deform[1]);
                if ( this->perform_timing_ ) { gt_timer_local_.stop(); }

                // ----------------------------------------
                // for every output N, warp the complex images
                // ----------------------------------------
                ArrayType complex_image(RO, E1, N, complex_image_raw.begin()+s*RO*E1*N);
                std::vector<ArrayType> warpped_complex_images_wider;

                if ( this->perform_timing_ ) { gt_timer_local_.start("perform_moco_warp_on_selected_images ... "); }
                this->perform_moco_warp_on_selected_images(selected_images_wider, complex_image, deform, warpped_complex_images_wider);
                if ( this->perform_timing_ ) { gt_timer_local_.stop(); }

                if ( !debug_folder_.empty() )
                {
                    for (n=0; n<dstN; n++)
                    {
                        std::stringstream os_local;
                        os_local << "_N_" << n << "_S_" << s;

                        gt_exporter_.export_array_complex(warpped_complex_images_wider[n], debug_folder_ + "warpped_complex_images" + os_local.str() + suffix_);
                    }
                }

                std::vector< std::vector<size_t> > loc_in_wider(dstN);
                for (n=0; n<dstN; n++)
                {
                    size_t num_of_images = selected_images[n].size();
                    size_t num_of_images_wider = selected_images_wider[n].size();

                    loc_in_wider[n].resize(num_of_images, 0);

                    size_t ind;
                    for (ii=0; ii<num_of_images; ii++)
                    {
                        size_t jj;
                        for (jj=0; jj<num_of_images_wider; jj++)
                        {
                            if(selected_images_wider[n][jj]==selected_images[n][ii])
                            {
                                loc_in_wider[n][ii] = jj;
                                break;
                            }
                        }
                    }
                }

                // ----------------------------------------
                // go back to multi-channel and fill the binning kspace
                // ----------------------------------------
                ArrayType kspace_binning(RO, E1, CHA, dstN, binning_obj_.kspace_binning_.begin()+s*RO*E1*CHA*dstN);
                ArrayType kspace_binning_wider(RO, E1, CHA, dstN, binning_obj_.kspace_binning_wider_.begin()+s*RO*E1*CHA*dstN);
                ArrayType kspace_binning_image_domain_average(RO, E1, CHA, dstN, binning_obj_.kspace_binning_image_domain_average_.begin()+s*RO*E1*CHA*dstN);
                hoNDArray< float > kspace_binning_hit_count(E1, dstN, binning_obj_.kspace_binning_hit_count_.begin()+s*E1*dstN);

                ArrayType coil_map(RO, E1, CHA, coil_map_raw.begin());
                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(coil_map, debug_folder_ + "coil_map" + os.str() + suffix_);

                ArrayType complexIm(RO, E1, CHA);

                ArrayType kspace_filled(RO, E1, CHA);
                hoNDArray< float > hit_count(E1);

                for (n=0; n<dstN; n++)
                {
                    GDEBUG_CONDITION_STREAM(this->verbose_, "Perform binning on step = " << n << " out of " << dstN);

                    ArrayType warpped_complex_images_multi_channel_wider;
                    ArrayType warpped_complex_images_multi_channel_image_domain_average_wider;

                    size_t num_selected_image_wider = selected_images_wider[n].size();

                    if(CHA>1)
                    {
                        // go back to multi-channel
                        if ( this->perform_timing_ ) { gt_timer_local_.start("go back to multi-channel ... "); }

                        warpped_complex_images_multi_channel_wider.create(RO, E1, CHA, num_selected_image_wider);

                        for (size_t ii=0; ii<num_selected_image_wider; ii++)
                        {
                            ArrayType complexIm2D(RO, E1, warpped_complex_images_wider[n].begin()+ii*RO*E1);
                            Gadgetron::multiply(coil_map, complexIm2D, complexIm);
                            memcpy(warpped_complex_images_multi_channel_wider.begin()+ii*RO*E1*CHA, complexIm.begin(), complexIm.get_number_of_bytes());
                        }

                        if ( this->perform_timing_ ) { gt_timer_local_.stop(); }
                    }
                    else
                    {
                        warpped_complex_images_multi_channel_wider.create(RO, E1, 1, num_selected_image_wider);

                        for (size_t ii=0; ii<num_selected_image_wider; ii++)
                        {
                            memcpy(warpped_complex_images_multi_channel_wider.begin()+ii*RO*E1, complex_image.begin() + selected_images_wider[n][ii]*RO*E1, complexIm.get_number_of_bytes());
                        }
                    }

                    if ( !debug_folder_.empty() )
                    {
                        std::stringstream os_local;
                        os_local << "_N_" << n << "_S_" << s;

                        gt_exporter_.export_array_complex(warpped_complex_images_multi_channel_wider, debug_folder_ + "warpped_complex_images_multi_channel" + os_local.str() + suffix_);
                    }

                    // average across all N
                    Gadgetron::sum_over_dimension(warpped_complex_images_multi_channel_wider, warpped_complex_images_multi_channel_image_domain_average_wider, 3);
                    Gadgetron::scal( (T)(1.0/num_selected_image_wider), warpped_complex_images_multi_channel_image_domain_average_wider);

                    if ( !debug_folder_.empty() )
                    {
                        std::stringstream os_local;
                        os_local << "_N_" << n << "_S_" << s;

                        gt_exporter_.export_array_complex(warpped_complex_images_multi_channel_image_domain_average_wider, debug_folder_ + "warpped_complex_images_multi_channel_image_domain_average" + os_local.str() + suffix_);
                    }

                    // go back to kspace
                    if ( this->perform_timing_ ) { gt_timer_local_.start("go back to kspace ... "); }
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(warpped_complex_images_multi_channel_wider);
                    if ( this->perform_timing_ ) { gt_timer_local_.stop(); }

                    // fill the binned kspace
                    if ( this->perform_timing_ ) { gt_timer_local_.start("fill the binned kspace ... "); }
                    this->fill_binned_kspace(s, n, selected_images_wider[n], warpped_complex_images_multi_channel_wider, kspace_filled, hit_count);
                    if ( this->perform_timing_ ) { gt_timer_local_.stop(); }

                    // copy results
                    memcpy(kspace_binning_wider.begin()+n*RO*E1*CHA, kspace_filled.begin(), kspace_filled.get_number_of_bytes());
                    memcpy(kspace_binning_image_domain_average.begin()+n*RO*E1*CHA, warpped_complex_images_multi_channel_image_domain_average_wider.begin(), warpped_complex_images_multi_channel_image_domain_average_wider.get_number_of_bytes());

                    // ---------------------------------------------

                    size_t num_of_images = selected_images[n].size();

                    ArrayType warpped_complex_images_multi_channel;
                    warpped_complex_images_multi_channel.create(RO, E1, CHA, num_of_images);

                    for (ii=0; ii<num_of_images; ii++)
                    {
                        memcpy(warpped_complex_images_multi_channel.begin()+ii*RO*E1*CHA, 
                            warpped_complex_images_multi_channel_wider.begin()+loc_in_wider[n][ii]*RO*E1*CHA, 
                            sizeof(std::complex<T>)*RO*E1*CHA);
                    }

                    this->fill_binned_kspace(s, n, selected_images[n], warpped_complex_images_multi_channel, kspace_filled, hit_count);

                    memcpy(kspace_binning.begin()+n*RO*E1*CHA, kspace_filled.begin(), kspace_filled.get_number_of_bytes());
                    memcpy(kspace_binning_hit_count.begin()+n*E1, hit_count.begin(), hit_count.get_number_of_bytes());

                    GDEBUG_CONDITION_STREAM(this->verbose_, "==================================================================");
                }

                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kspace_binning_image_domain_average, debug_folder_ + "kspace_binning_image_domain_average_IMAGE" + os.str() + suffix_);

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(kspace_binning_image_domain_average);

                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kspace_binning, debug_folder_ + "kspace_binning" + os.str() + suffix_);
                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kspace_binning_wider, debug_folder_ + "kspace_binning_wider" + os.str() + suffix_);
                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kspace_binning_image_domain_average, debug_folder_ + "kspace_binning_image_domain_average" + os.str() + suffix_);
                if ( !debug_folder_.empty() ) gt_exporter_.export_array(kspace_binning_hit_count, debug_folder_ + "kspace_binning_hit_count" + os.str() + suffix_);
            }
        }
        else
        {
            // to be implemented
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::compute_kspace_binning() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::interpolate_best_HB_images(const std::vector<float>& cpt_time_ratio_bestHB, const hoNDArray<T>& mag_bestHB, const std::vector<float>& desired_cpt, hoNDArray<T>& mag_bestHB_at_desired_cpt)
{
    try
    {
        size_t RO = mag_bestHB.get_size(0);
        size_t E1 = mag_bestHB.get_size(1);
        size_t N = mag_bestHB.get_size(2);

        size_t dstN = desired_cpt.size();

        mag_bestHB_at_desired_cpt.create( RO, E1, dstN );

        // check the boundary of cpt_time_ratio_bestHB
        std::vector<float> cpt_time_ratio_bestHB_checked(cpt_time_ratio_bestHB);

        if(cpt_time_ratio_bestHB_checked[0] > cpt_time_ratio_bestHB_checked[1])
        {
            cpt_time_ratio_bestHB_checked[0] = 1.0 - cpt_time_ratio_bestHB_checked[0];
            if(cpt_time_ratio_bestHB_checked[0] > cpt_time_ratio_bestHB_checked[1])
            {
                cpt_time_ratio_bestHB_checked[0] = 0;
            }

            GDEBUG_STREAM("Cpt time ratio 0 used for interpolating : " << cpt_time_ratio_bestHB_checked[0]);
        }

        if(cpt_time_ratio_bestHB_checked[N-1] < cpt_time_ratio_bestHB_checked[N-2])
        {
            cpt_time_ratio_bestHB_checked[N-1] = 1.0;
            GDEBUG_STREAM("Cpt time ratio N-1 used for interpolating : " << cpt_time_ratio_bestHB_checked[N-1]);
        }

        // extend the best HB for periodic boundary condition
        hoNDArray<T> mag(RO, E1, N+2);
        std::vector<float> cpt(N+2);

        memcpy(mag.begin()+RO*E1, mag_bestHB.begin(), mag_bestHB.get_number_of_bytes());
        std::copy(cpt_time_ratio_bestHB_checked.begin(), cpt_time_ratio_bestHB_checked.end(), cpt.begin() + 1);

        size_t n;

        // first phase
        memcpy(mag.begin(), mag_bestHB.begin(), sizeof(T)*RO*E1);
        cpt[0] = (float)( -(1.0 - cpt_time_ratio_bestHB_checked[N-1]) );

        // last phase
        // memcpy(mag.begin()+(N+1)*RO*E1, mag_bestHB.begin(), sizeof(T)*RO*E1);
        memcpy(mag.begin()+(N+1)*RO*E1, mag_bestHB.begin()+(N-1)*RO*E1, sizeof(T)*RO*E1);
        cpt[N+1] = (float)( 1.0 + cpt_time_ratio_bestHB_checked[0] );

        if ( !debug_folder_.empty() ) gt_exporter_.export_array(mag, debug_folder_ + "mag_bestHB_periodic_condition" + suffix_);

        GDEBUG_CONDITION_STREAM(this->verbose_, "Cpt time ratio used for interpolating best heart beat : ");
        for ( n=0; n<N+2; n++ )
        {
            GDEBUG_CONDITION_STREAM(this->verbose_, " n - " << (long long)(n)-1 << " - cpt time ratio : " << cpt[n]);
        }

        GDEBUG_CONDITION_STREAM(this->verbose_, "Interpolate the best cardiac cycle using bspline ... ");

        size_t SplineDegree = 5;

        long long ro, e1;

        std::vector<double> x(N+2);
        for ( n=0; n<N+2; n++ ) x[n] = cpt[n];

#pragma omp parallel default(none) shared(x, SplineDegree, RO, E1, N, mag, dstN, desired_cpt, mag_bestHB_at_desired_cpt) private(ro, e1, n)
        {
            hoNDArray<T> y(N+2);
            hoNDArray<T> coeff(N+2);
            hoNDBSpline<T, 1 > interp;

            for ( e1=0; e1<E1; e1++ )
            {
                for ( ro=0; ro<RO; ro++ )
                {
                    for ( n=0; n<N+2; n++ )
                    {
                        y(n) = mag(ro, e1, n);
                    }

                    // compute the coefficient
                    interp.computeBSplineCoefficients(y, SplineDegree, coeff);

                    for ( n=0; n<dstN; n++ )
                    {
                        float dst_n = (N+1)*(desired_cpt[n]-x[0])/(x[N+1] - x[0]);
                        mag_bestHB_at_desired_cpt(ro, e1, n) = interp.evaluateBSpline(coeff.begin(), N+2, SplineDegree, 0, dst_n);
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::interpolate_best_HB_images() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::select_images_with_navigator(float desired_cpt_time_ratio, float accepted_nav[2], size_t n, size_t s, std::vector<size_t>& selected)
{
    try
    {
        size_t RO = binning_obj_.complex_image_raw_.get_size(0);
        size_t E1 = binning_obj_.complex_image_raw_.get_size(1);

        size_t N = binning_obj_.complex_image_raw_.get_size(3);

        // find all images including kspace lines falling into the allowed temporal ratio window
        float cpt_ratio_step_size = binning_obj_.desired_cpt_[1] - binning_obj_.desired_cpt_[0];

        float accepted_cpt[2];
        accepted_cpt[0] = desired_cpt_time_ratio - this->kspace_binning_max_temporal_window_*cpt_ratio_step_size;
        accepted_cpt[1] = desired_cpt_time_ratio + this->kspace_binning_max_temporal_window_*cpt_ratio_step_size;

        float accepted_cpt_mirror[2];
        accepted_cpt_mirror[0] = accepted_cpt[0];
        accepted_cpt_mirror[1] = accepted_cpt[1];

        if ( accepted_cpt[0] < 0 )
        {
            float delta = -accepted_cpt[0];
            accepted_cpt[0] = 0.0f;

            accepted_cpt_mirror[0] = 1.0f -delta;
            accepted_cpt_mirror[1] = 1.0f;
        }

        if ( accepted_cpt[1] > 1.0 )
        {
            float delta = (float)(accepted_cpt[1] - 1.0);
            accepted_cpt[1] = 1.0f;

            accepted_cpt_mirror[0] = 0.0f;
            accepted_cpt_mirror[1] = delta;
        }

        float temporal_window = this->kspace_binning_max_temporal_window_;

        size_t e1, nn;
        while ( (selected.size() < binning_obj_.accel_factor_E1_) && (temporal_window<6.0) )
        {
            selected.clear();

            for ( nn=0; nn<N; nn++ )
            {
                bool cpt_in_range = false;
                for (e1=0; e1<E1; e1++)
                {
                    float cpt = binning_obj_.cpt_time_ratio_(e1, nn, s);
                    if(cpt>0)
                    {
                        if ( (cpt>=accepted_cpt[0] && cpt<accepted_cpt[1]) || (cpt>=accepted_cpt_mirror[0] && cpt<accepted_cpt_mirror[1]) )
                        {
                            cpt_in_range = true;
                            break;
                        }
                    }
                }

                if(cpt_in_range)
                {
                    float nav = binning_obj_.navigator_(E1/2, nn, s);
                    if ( nav>=accepted_nav[0] && nav<=accepted_nav[1] )
                    {
                        selected.push_back(nn);
                    }
                }
            }

            temporal_window *= 1.1f;

            accepted_cpt[0] = desired_cpt_time_ratio - this->kspace_binning_max_temporal_window_/2*cpt_ratio_step_size;
            accepted_cpt[1] = desired_cpt_time_ratio + this->kspace_binning_max_temporal_window_/2*cpt_ratio_step_size;
            accepted_cpt_mirror[0] = accepted_cpt[0];
            accepted_cpt_mirror[1] = accepted_cpt[1];

            if ( accepted_cpt[0] < 0 )
            {
                float delta = -accepted_cpt[0];
                accepted_cpt[0] = 0.0f;

                accepted_cpt_mirror[0] = 1.0f -delta;
                accepted_cpt_mirror[1] = 1.0f;
            }

            if ( accepted_cpt[1] > 1.0 )
            {
                float delta = (float)(accepted_cpt[1] - 1.0);
                accepted_cpt[1] = 1.0f;

                accepted_cpt_mirror[0] = 0.0f;
                accepted_cpt_mirror[1] = delta;
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::select_images_with_navigator() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::perform_moco_selected_images_with_best_heart_beat(const std::vector < std::vector<size_t> >& selected_images, const hoNDArray<T>& mag_s, const hoNDArray<T>& mag_bestHB_at_desired_cpt, DeformationFieldContinerType& dx, DeformationFieldContinerType& dy)
{
    try
    {
        size_t RO = mag_s.get_size(0);
        size_t E1 = mag_s.get_size(1);
        size_t N = mag_s.get_size(2);

        size_t dstN = selected_images.size();

        size_t n, p;

        ImageContinerType input;
        std::vector<size_t> cols(dstN);

        for ( n=0; n<dstN; n++ )
        {
            cols[n] = selected_images[n].size() + 1;
        }

        std::vector<size_t> dim(2);
        dim[0] = RO;
        dim[1] = E1;

        input.create(cols, true);

        for ( n=0; n<dstN; n++ )
        {
            size_t num = selected_images[n].size();

            input(n, 0).create(dim, const_cast<T*>(mag_bestHB_at_desired_cpt.begin()+n*RO*E1), false);

            for ( p=0; p<num; p++ )
            {
                input(n, p+1).create(dim, const_cast<T*>(mag_s.begin()+selected_images[n][p]*RO*E1), false);
            }

            if ( !debug_folder_.empty() )
            {
                std::ostringstream ostr;
                ostr << "Picked_images_before_moco_n" << n;
                std::string fileNameUsed = ostr.str();

                hoNDArray<T> tmp;
                input.to_NDArray(n, tmp);
                gt_exporter_.export_array(tmp, debug_folder_+fileNameUsed + suffix_);
            }
        }

        std::vector<unsigned int> key_frame(dstN, 0);

        Gadgetron::hoImageRegContainer2DRegistration<hoNDImage<float, 2>, hoNDImage<float, 2>, double> reg;

        GDEBUG_STREAM("Perform moco against best heart beat : " << this->kspace_binning_moco_reg_strength_);
        GDEBUG_STREAM("MOCO iterations : ");
        for (size_t ii=0; ii<this->kspace_binning_moco_iters_.size(); ii++)
        {
            GDEBUG_STREAM(this->kspace_binning_moco_iters_[ii]);
        }

        bool warp_image = false;
        if ( !debug_folder_.empty() ) warp_image = true;

        Gadgetron::perform_moco_fixed_key_frame_2DT(input, key_frame, 
            this->kspace_binning_moco_reg_strength_, this->kspace_binning_moco_iters_, false, warp_image, reg);

        dx.copyFrom(reg.deformation_field_[0]);
        dy.copyFrom(reg.deformation_field_[1]);

        if(warp_image && !debug_folder_.empty())
        {
            hoNDArray<T> im;

            for ( n=0; n<dstN; n++ )
            {
                reg.warped_container_.to_NDArray(n, im);
                {
                    std::ostringstream ostr;
                    ostr << "Picked_images_after_moco_n" << n;
                    std::string fileNameUsed = ostr.str();

                    gt_exporter_.export_array(im, debug_folder_+fileNameUsed + suffix_);
                }
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::perform_moco_selected_images_with_best_heart_beat() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::perform_moco_warp_on_selected_images( const std::vector < std::vector<size_t> >& selected_images, const ArrayType& complex_image, DeformationFieldContinerType deform[2], std::vector<ArrayType>& warpped_complex_images)
{
    try
    {
        size_t RO = complex_image.get_size(0);
        size_t E1 = complex_image.get_size(1);
        size_t N = complex_image.get_size(2);

        size_t dstN = selected_images.size();

        size_t n, p;

        std::vector<size_t> cols(dstN);

        for ( n=0; n<dstN; n++ )
        {
            cols[n] = selected_images[n].size() + 1;
        }

        std::vector<size_t> dim(2);
        dim[0] = RO;
        dim[1] = E1;

        ComplexImageContinerType input;
        input.create(cols, true);

        for ( n=0; n<dstN; n++ )
        {
            size_t num = selected_images[n].size();

            input(n, 0).create(dim);
            Gadgetron::clear(input(n, 0));

            for ( p=0; p<num; p++ )
            {
                input(n, p+1).create(dim, const_cast< std::complex<T>* >(complex_image.begin()+selected_images[n][p]*RO*E1), false);
            }
        }

        ComplexImageContinerType output;
        output.copyFrom(input);

        Gadgetron::hoImageRegContainer2DRegistration<ImageType, ImageType, double> reg;
        reg.warpContainer2D(input, input, deform, output);

        warpped_complex_images.resize(dstN);
        for ( n=0; n<dstN; n++ )
        {
            size_t num = selected_images[n].size();
            warpped_complex_images[n].create(RO, E1, num);

            for (size_t ii=0; ii<num; ii++)
            {
                memcpy(warpped_complex_images[n].begin()+ii*RO*E1, output(n, ii+1).begin(), output(n, ii+1).get_number_of_bytes());
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::perform_moco_warp_on_selected_images() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::fill_binned_kspace(size_t s, size_t dst_n, const std::vector<size_t>& selected_images, const ArrayType& warpped_kspace, ArrayType& kspace_filled, hoNDArray<float>& hit_count)
{
    try
    {
        size_t RO = warpped_kspace.get_size(0);
        size_t E1 = warpped_kspace.get_size(1);
        size_t CHA = warpped_kspace.get_size(2);
        size_t N = warpped_kspace.get_size(3);

        kspace_filled.create(RO, E1, CHA);
        Gadgetron::clear(kspace_filled);

        hit_count.create(E1);
        Gadgetron::clear(hit_count);

        float desired_cpt_time_ratio = binning_obj_.desired_cpt_[dst_n];

        float cpt_time_ratio_step_size = binning_obj_.desired_cpt_[1] - binning_obj_.desired_cpt_[0];

        float accepted_cpt_time_ratio[2], accepted_cpt_time_ratio_mirror[2];

        long long required_E1_range[2];
        long long required_E1_block = 1;
        if ( required_E1_block > 0 )
        {
            required_E1_range[0] = (long long)(E1/2 - required_E1_block*binning_obj_.accel_factor_E1_/2);
            required_E1_range[1] = (long long)(required_E1_range[0] + required_E1_block*binning_obj_.accel_factor_E1_/2 - 1);
        }
        else
        {
            required_E1_range[0] = 0;
            required_E1_range[1] = (long long)(E1-1);
        }

        long long e1, cha, n;

        bool central_lines_filled = false;
        long long selectedN = static_cast<long long>(selected_images.size());

        // if the picked phases are not enough, enlarge the temporal window
        ArrayType readOut1(RO, CHA);
        ArrayType readOut2(RO, CHA);
        ArrayType readOut3(RO, CHA);

        // if the minimal cardiac phase temporal window is set, adjust the accepted_cpt_time_ratio range
        float minimal_cpt_time_ratio = -1;
        if ( this->kspace_binning_minimal_cardiac_phase_width_ > 0 )
        {
            minimal_cpt_time_ratio = this->kspace_binning_minimal_cardiac_phase_width_/binning_obj_.mean_RR_;
        }

        float cpt_time_ratio_window = 0.5;
        float meanRR = binning_obj_.mean_RR_;

        while ( !central_lines_filled && (cpt_time_ratio_window<this->kspace_binning_max_temporal_window_) )
        {
            accepted_cpt_time_ratio[0] = desired_cpt_time_ratio - cpt_time_ratio_window*cpt_time_ratio_step_size;
            accepted_cpt_time_ratio[1] = desired_cpt_time_ratio + cpt_time_ratio_window*cpt_time_ratio_step_size;

            if ( minimal_cpt_time_ratio > 0 )
            {
                float currentRange = accepted_cpt_time_ratio[1]-accepted_cpt_time_ratio[0];

                if ( currentRange < minimal_cpt_time_ratio )
                {
                    // enlarge the accepted_cpt_time_ratio range
                    float diff = (minimal_cpt_time_ratio-currentRange)/2;

                    accepted_cpt_time_ratio[0] -= diff;
                    accepted_cpt_time_ratio[1] += diff;

                    GDEBUG_CONDITION_STREAM(true, "accepted cardiac phase ratio is increased to " 
                        << minimal_cpt_time_ratio << " - meanRR - " << meanRR << " - minimal temporal window - " 
                        << this->kspace_binning_minimal_cardiac_phase_width_/1000);
                }
            }

            accepted_cpt_time_ratio_mirror[0] = accepted_cpt_time_ratio[0];
            accepted_cpt_time_ratio_mirror[1] = accepted_cpt_time_ratio[1];

            if ( accepted_cpt_time_ratio[0] < 0 )
            {
                float delta = -accepted_cpt_time_ratio[0];
                accepted_cpt_time_ratio[0] = 0.0f;

                accepted_cpt_time_ratio_mirror[0] = 1.0f -delta;
                accepted_cpt_time_ratio_mirror[1] = 1.0f;
            }

            if ( accepted_cpt_time_ratio[1] > 1.0 )
            {
                float delta = accepted_cpt_time_ratio[1] - 1.0f;
                accepted_cpt_time_ratio[1] = 1.0f;

                accepted_cpt_time_ratio_mirror[0] = 0.0f;
                accepted_cpt_time_ratio_mirror[1] = delta;
            }

            GDEBUG_CONDITION_STREAM(this->verbose_, "Fill kspace - accepted_cpt_time_ratio : " << accepted_cpt_time_ratio[0] << " - " << accepted_cpt_time_ratio[1]);
            GDEBUG_CONDITION_STREAM(this->verbose_, "Fill kspace - accepted_cpt_time_ratio_mirror : " << accepted_cpt_time_ratio_mirror[0] << " - " << accepted_cpt_time_ratio_mirror[1]);

            for ( e1=0; e1<E1; e1++ )
            {
                hit_count(e1) = 0.0f;
            }

            Gadgetron::clear(kspace_filled);

            std::vector<long long> lineFilled;

            for ( n=0; n<selectedN; n++ )
            {
                for ( e1=0; e1<E1; e1++ )
                {
                    float cptLIN = binning_obj_.cpt_time_ratio_(e1, selected_images[n], s);

                    if ( (cptLIN>=accepted_cpt_time_ratio[0] && cptLIN<=accepted_cpt_time_ratio[1]) || (cptLIN>=accepted_cpt_time_ratio_mirror[0] && cptLIN<=accepted_cpt_time_ratio_mirror[1]) )
                    {
                        for ( cha=0; cha<CHA; cha++ )
                        {
                            memcpy(readOut1.begin()+cha*RO, kspace_filled.begin()+cha*RO*E1+e1*RO, sizeof(std::complex<T>)*RO);

                            memcpy(readOut2.begin()+cha*RO, warpped_kspace.begin()+n*RO*E1*CHA+cha*RO*E1+e1*RO, sizeof(std::complex<T>)*RO);
                        }

                        Gadgetron::add(readOut1, readOut2, readOut3);

                        for ( cha=0; cha<CHA; cha++ )
                        {
                            memcpy(kspace_filled.begin()+cha*RO*E1+e1*RO, readOut3.begin()+cha*RO, sizeof(std::complex<T>)*RO);
                        }

                        // increase the counter
                        float value = hit_count(e1);
                        hit_count(e1) = value+1;

                        // check whether the central line is filled
                        if ( e1>=required_E1_range[0] && e1<=required_E1_range[1] )
                        {
                            central_lines_filled = true;
                        }

                        lineFilled.push_back(e1);
                    }
                }
            }

            // print out the line filled
            GDEBUG_CONDITION_STREAM(this->verbose_, "Destination " << dst_n << " - temporal window ratio " << cpt_time_ratio_window << " - number of lines filled : " << lineFilled.size());

            for ( e1=0; e1<E1; e1++ )
            {
                if ( hit_count(e1) > 1.0f )
                {
                    for ( cha=0; cha<CHA; cha++ )
                    {
                        memcpy(readOut1.begin()+cha*RO, 
                            kspace_filled.begin()+cha*RO*E1+e1*RO,
                            sizeof(std::complex<T>)*RO);
                    }

                    Gadgetron::scal(1.0f/hit_count(e1), readOut1);

                    for ( cha=0; cha<CHA; cha++ )
                    {
                        memcpy(kspace_filled.begin()+cha*RO*E1+e1*RO, 
                            readOut1.begin()+cha*RO, sizeof(std::complex<T>)*RO);
                    }
                }
            }

            if ( !central_lines_filled ) cpt_time_ratio_window *= 1.1f;
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::fill_binned_kspace() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::perform_recon_binned_kspace(const std::vector<size_t>& slices_not_processing)
{
    try
    {
        hoNDArray< float >& kspace_binning_hit_count = binning_obj_.kspace_binning_hit_count_;
        ArrayType& coil_map = binning_obj_.coil_map_raw_;

        ArrayType& kspace_binning = binning_obj_.kspace_binning_;
        ArrayType& kspace_binning_wider = binning_obj_.kspace_binning_wider_;
        ArrayType& kspace_binning_image_domain_average = binning_obj_.kspace_binning_image_domain_average_;

        ArrayType& complex_image_binning = binning_obj_.complex_image_binning_;

        size_t RO = kspace_binning.get_size(0);
        size_t E1 = kspace_binning.get_size(1);
        size_t CHA = kspace_binning.get_size(2);
        size_t N = kspace_binning.get_size(3);
        size_t S = kspace_binning.get_size(4);

        complex_image_binning.create(RO, E1, 1, N, S);
        Gadgetron::clear(complex_image_binning);

        size_t e1, cha, n, s;

        if(CHA==1)
        {
            // use the neighboring frames if there are holes in the binned kspace
            for (s=0; s<S; s++)
            {
                bool not_processing = false;
                for (size_t kk=0; kk<slices_not_processing.size(); kk++)
                {
                    if(slices_not_processing[kk] == s)
                    {
                        not_processing = true;
                        break;
                    }
                }

                if(not_processing)
                {
                    GWARN_STREAM("Due to previously happened errors, slice " << s << " will not be processed ... ");
                    continue;
                }

                for (n=0; n<N; n++)
                {
                    for (e1=0; e1<E1; e1++)
                    {
                        if(kspace_binning_hit_count(e1, n, s)==0)
                        {
                            long long ind_left(0), n_left(n), ind_right(0), n_right(n);

                            ind_left = 0;
                            while (kspace_binning_hit_count(e1, n_left, s)==0 && (ind_left<N))
                            {
                                n_left--;
                                if(n_left<0) n_left += N;
                                ind_left++;
                            }

                            ind_right = 0;
                            while (kspace_binning_hit_count(e1, n_right, s)==0 && (ind_right<N))
                            {
                                n_right++;
                                if(n_right>=N) n_right -= N;
                                ind_right++;
                            }

                            if(ind_left>=N && ind_right>=N)
                            {
                                // cannot fill the hole ...
                            }
                            else if(ind_left<N && ind_right>=N)
                            {
                                // use the left
                                memcpy(&kspace_binning(0, e1, n, s), &kspace_binning(0, e1, n_left, s), sizeof(std::complex<T>)*RO);
                            }
                            else if(ind_left>=N && ind_right<N)
                            {
                                // use the right
                                memcpy(&kspace_binning(0, e1, n, s), &kspace_binning(0, e1, n_right, s), sizeof(std::complex<T>)*RO);
                            }
                            else if(ind_left<N && ind_right<N)
                            {
                                ArrayType readout_left(RO);
                                memcpy(readout_left.begin(), &kspace_binning(0, e1, n_left, s), sizeof(std::complex<T>)*RO);

                                ArrayType readout_right(RO);
                                memcpy(readout_right.begin(), &kspace_binning(0, e1, n_right, s), sizeof(std::complex<T>)*RO);

                                T w_l = (T)(ind_right)/(ind_left+ind_right);
                                T w_r = (T)(ind_left)/(ind_left+ind_right);

                                Gadgetron::scal(w_l, readout_left);
                                Gadgetron::scal(w_r, readout_right);

                                Gadgetron::add(readout_left, readout_right, readout_left);

                                memcpy(&kspace_binning(0, e1, n, s), readout_left.begin(), sizeof(std::complex<T>)*RO);
                            }
                        }
                    }
                }
            }

            if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kspace_binning, debug_folder_ + "kspace_binning_before_fft_recon" + suffix_);

            // perform fft recon
            ArrayType reconFFT(kspace_binning);
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspace_binning, complex_image_binning);
        }
        else
        {
            for (s=0; s<S; s++)
            {
                bool not_processing = false;
                for (size_t kk=0; kk<slices_not_processing.size(); kk++)
                {
                    if(slices_not_processing[kk] == s)
                    {
                        not_processing = true;
                        break;
                    }
                }

                if(not_processing)
                {
                    GWARN_STREAM("Due to previously happened errors, slice " << s << " will not be processed ... ");
                    continue;
                }

                std::stringstream os;
                os << "_S_" << s;

                ArrayType kspace(RO, E1, CHA, N, 1, kspace_binning.begin()+s*RO*E1*CHA*N);
                ArrayType kspace_wider(RO, E1, CHA, N, 1, kspace_binning_wider.begin()+s*RO*E1*CHA*N);
                ArrayType kspaceRef(RO, E1, CHA, N, 1, kspace_binning_image_domain_average.begin()+s*RO*E1*CHA*N);
                ArrayType coilMap(RO, E1, CHA, coil_map.begin());

                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kspace, debug_folder_ + "kspace_binning_linear_recon_kspace" + os.str() + suffix_);
                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kspace_wider, debug_folder_ + "kspace_binning_linear_recon_kspace_wider" + os.str() + suffix_);
                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kspaceRef, debug_folder_ + "kspace_binning_linear_recon_kspaceRef" + os.str() + suffix_);
                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(coilMap, debug_folder_ + "kspace_binning_linear_recon_coilMap" + os.str() + suffix_);

                // perform linear recon
                ArrayType resKSpace, resIm, kernel, kernelIm;

                if ( this->perform_timing_ ) { gt_timer_local_.start("perform linear recon on kspace binning ... "); }

                if(this->use_nonlinear_binning_recon_)
                {
                    this->perform_linear_recon_on_kspace_binning(kspace_wider, kspaceRef, coilMap, resKSpace, resIm, kernel, kernelIm);
                }
                else
                {
                    this->perform_linear_recon_on_kspace_binning(kspace, kspaceRef, coilMap, resKSpace, resIm, kernel, kernelIm);
                }
                if ( this->perform_timing_ ) { gt_timer_local_.stop(); }

                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(resKSpace, debug_folder_ + "kspace_binning_linear_recon_resKSpace" + os.str() + suffix_);
                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(resIm, debug_folder_ + "kspace_binning_linear_recon_resIm" + os.str() + suffix_);
                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kernel, debug_folder_ + "kspace_binning_linear_recon_kernel" + os.str() + suffix_);
                if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kernelIm, debug_folder_ + "kspace_binning_linear_recon_kernelIm" + os.str() + suffix_);

                // perform nonlinear recon
                if(this->use_nonlinear_binning_recon_)
                {
                    ArrayType resKSpaceNonLinear, resImNonLinear;

                    if ( this->perform_timing_ ) { gt_timer_local_.start("perform non-linear recon on kspace binning ... "); }
                    this->perform_non_linear_recon_on_kspace_binning(kspace, resKSpace, coilMap, kernel, kernelIm, resKSpaceNonLinear, resImNonLinear);
                    if ( this->perform_timing_ ) { gt_timer_local_.stop(); }

                    if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(resKSpaceNonLinear, debug_folder_ + "kspace_binning_linear_recon_resKSpaceNonLinear" + os.str() + suffix_);
                    if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(resImNonLinear, debug_folder_ + "kspace_binning_linear_recon_resImNonLinear" + os.str() + suffix_);

                    memcpy(complex_image_binning.begin()+s*RO*E1*N, resImNonLinear.begin(), resImNonLinear.get_number_of_bytes());
                }
                else
                {
                    memcpy(complex_image_binning.begin()+s*RO*E1*N, resIm.begin(), resIm.get_number_of_bytes());
                }
            }
        }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::perform_recon_binned_kspace() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::perform_linear_recon_on_kspace_binning(const ArrayType& kspace, const ArrayType& kspaceInitial, const ArrayType& coilMap, ArrayType& resKSpace, ArrayType& resIm, ArrayType& kernel, ArrayType& kernelIm)
{
    try
    {
        Gadgetron::GadgetronTimer timer;
        timer.set_timing_in_destruction(false);

        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t CHA = kspace.get_size(2);
        size_t N = kspace.get_size(3);

        size_t startE1, endE1;
        startE1 = this->binning_obj_.sampling_.sampling_limits_[1].min_;
        endE1 = this->binning_obj_.sampling_.sampling_limits_[1].max_;

        // average all N
        ArrayType acs;
        Gadgetron::sum_over_dimension(kspaceInitial, acs, 3);

        // estimate kernel
        size_t kRO = kspace_binning_kSize_RO_;
        size_t kE1 = kspace_binning_kSize_E1_;

        size_t convKRO = 2 * kRO - 1;
        size_t convKE1 = 2 * kE1 - 1;

        double reg_lamda = kspace_binning_reg_lamda_;
        double over_determine_ratio = 45;

        ArrayType acsSrc(RO, E1, CHA, acs.begin());
        ArrayType acsDst(RO, E1, CHA, acs.begin());

        kernel.create(convKRO, convKE1, CHA, CHA);
        kernelIm.create(RO, E1, CHA, CHA);

        if ( this->perform_timing_ ) { timer.start("spirit calibration ... "); }
        Gadgetron::spirit2d_calib_convolution_kernel(acsSrc, acsDst, reg_lamda, kRO, kE1, 1, 1, kernel, true);
        if ( this->perform_timing_ ) { timer.stop(); }

        if ( this->perform_timing_ ) { timer.start("spirit image domain kernel ... "); }
        Gadgetron::spirit2d_image_domain_kernel(kernel, RO, E1, kernelIm);
        if ( this->perform_timing_ ) { timer.stop(); }

        // perform recon
        size_t iter_max = kspace_binning_linear_iter_max_;
        double iter_thres = kspace_binning_linear_iter_thres_;
        bool print_iter = true;

        ArrayType emptyInitial;
        if ( this->perform_timing_ ) { timer.start("spirit linear recon ... "); }
        Gadgetron::perform_spirit_recon_linear_2DT(kspace, startE1, endE1, kernelIm, emptyInitial, resKSpace, iter_max, iter_thres, print_iter);
        if ( this->perform_timing_ ) { timer.stop(); }

        if ( this->perform_timing_ ) { timer.start("coil combination ... "); }
        ArrayType complexIm;
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(resKSpace, complexIm);
        Gadgetron::coil_combine(complexIm, coilMap, 2, resIm);
        if ( this->perform_timing_ ) { timer.stop(); }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::perform_linear_recon_on_kspace_binning() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::perform_non_linear_recon_on_kspace_binning(const ArrayType& kspace, const ArrayType& kspaceLinear, const ArrayType& coilMap, const ArrayType& kernel, const ArrayType& kernelIm, ArrayType& resKSpace, ArrayType& resIm)
{
    try
    {
        Gadgetron::GadgetronTimer timer;
        timer.set_timing_in_destruction(false);

        bool print_iter = true;

        if ( this->perform_timing_ ) { timer.start("spirit non linear recon ... "); }

        Gadgetron::perform_spirit_recon_non_linear_2DT(kspace, kernelIm, coilMap, kspaceLinear, resKSpace, 
                                                    kspace_binning_nonlinear_iter_max_, 
                                                    kspace_binning_nonlinear_iter_thres_, 
                                                    kspace_binning_nonlinear_data_fidelity_lamda_, 
                                                    kspace_binning_nonlinear_image_reg_lamda_, 
                                                    kspace_binning_nonlinear_reg_N_weighting_ratio_, 
                                                    kspace_binning_nonlinear_reg_use_coil_sen_map_, 
                                                    kspace_binning_nonlinear_reg_with_approx_coeff_, 
                                                    kspace_binning_nonlinear_reg_wav_name_, 
                                                    print_iter);

        if ( this->perform_timing_ ) { timer.stop(); }

        if ( this->perform_timing_ ) { timer.start("coil combination ... "); }
        ArrayType complexIm;
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(resKSpace, complexIm);
        Gadgetron::coil_combine(complexIm, coilMap, 2, resIm);
        if ( this->perform_timing_ ) { timer.stop(); }
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::perform_non_linear_recon_on_kspace_binning() ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCMR CmrKSpaceBinning< float >;

}
