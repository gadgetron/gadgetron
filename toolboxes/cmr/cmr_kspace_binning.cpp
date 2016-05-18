/** \file   cmr_kspace_binning.cpp
    \brief  Implement kspace binning recon for 2D acquisition
    \author Hui Xue
*/

#include "cmr_kspace_binning.h"
#include "log.h"
#include "hoNDFFT.h"
#include "mri_core_grappa.h"
#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_coil_map_estimation.h"

#include "cmr_time_stamp.h"
#include "cmr_motion_correction.h"

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

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kspace, debug_folder_ + "binning_kspace");

        // -----------------------------------------------------
        // perform the raw data recon
        // -----------------------------------------------------
        if ( this->perform_timing_ ) { gt_timer_.start("raw data recon ... "); }

        this->perform_raw_data_recon();

        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(binning_obj_.full_kspace_raw_, debug_folder_ + "full_kspace_raw");
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(binning_obj_.complex_image_raw_, debug_folder_ + "complex_image_raw");
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(binning_obj_.coil_map_raw_, debug_folder_ + "coil_map_raw");

        // -----------------------------------------------------
        // estimate acquisition time stamp and cardiac phase time ratio for all acuqired kspace lines
        // count the heart beats and when it starts
        // -----------------------------------------------------
        binning_obj_.cpt_time_stamp_.create(E1, N, S);
        binning_obj_.ind_heart_beat_.create(E1, N, S);

        if ( this->perform_timing_ ) { gt_timer_.start("estimate time stamps ... "); }

        this->estimate_time_stamps();

        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.time_stamp_, debug_folder_ + "time_stamp");
        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.cpt_time_stamp_, debug_folder_ + "cpt_time_stamp");
        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.phs_time_stamp_, debug_folder_ + "phs_time_stamp");
        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.phs_cpt_time_stamp_, debug_folder_ + "phs_cpt_time_stamp");

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

        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.navigator_, debug_folder_ + "respiratory_navigator");

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
        // perform kspace binning
        // -----------------------------------------------------
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
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(data, debug_folder_ + "raw_data_recon_data");

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

        Gadgetron::compute_averaged_data_N_S(data, average_all_ref_N, average_all_ref_S, count_sampling_freq, ref);

        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(ref, debug_folder_ + "raw_data_recon_ref");

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
                T v = Gadgetron::norm2(dataCha);
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

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(ref_src, debug_folder_ + "raw_data_recon_ref_src");
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(ref_dst, debug_folder_ + "raw_data_recon_ref_dst");

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

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(filter_RO_ref_coi_map, debug_folder_ + "raw_data_recon_filter_RO_ref_coi_map");
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(filter_E1_ref_coi_map, debug_folder_ + "raw_data_recon_filter_E1_ref_coi_map");

        // apply ref coil map filter
        ArrayType ref_coil_map_dst;
        Gadgetron::apply_kspace_filter_ROE1(ref_dst, filter_RO_ref_coi_map, filter_E1_ref_coi_map, ref_coil_map_dst);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(ref_coil_map_dst, debug_folder_ + "raw_data_recon_ref_coil_map_dst");

        // compute coil map
        ArrayType complex_im_coil_map;
        Gadgetron::hoNDFFT<T>::instance()->ifft2c(ref_coil_map_dst, complex_im_coil_map);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(complex_im_coil_map, debug_folder_ + "raw_data_recon_complex_im_coil_map");

        size_t ks = 7;
        size_t power = 3;

        if ( this->perform_timing_ ) { gt_timer_.start("--> perform_raw_data_recon, compute coil map"); }

        Gadgetron::coil_map_2d_Inati(complex_im_coil_map, coil_map, ks, power);

        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(coil_map, debug_folder_ + "raw_data_recon_Result_coil_map");

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

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(convKer, debug_folder_ + "raw_data_recon_convKer");
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(kernelIm, debug_folder_ + "raw_data_recon_kernelIm");

        // compute unmixing coefficient
        ArrayType unmixing_coeff;
        hoNDArray<T> gFactor;
        Gadgetron::grappa2d_unmixing_coeff(kernelIm, coil_map, (size_t)binning_obj_.accel_factor_E1_, unmixing_coeff, gFactor);

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(unmixing_coeff, debug_folder_ + "raw_data_recon_unmixing_coeff");
        if ( !debug_folder_.empty() ) gt_exporter_.export_array(gFactor, debug_folder_ + "raw_data_recon_gFactor");

        // ------------------------------------------------
        // perform reconstruction
        // ------------------------------------------------

        // compute aliased images
        ArrayType aliased_im(data);
        Gadgetron::hoNDFFT<T>::instance()->ifft2c(data, aliased_im);
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(aliased_im, debug_folder_ + "raw_data_recon_aliased_im");

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
        }
        else
        {
            // compute single-channel full kspace
            Gadgetron::hoNDFFT<T>::instance()->fft2c(complex_image, full_kspace);
        }

        if ( this->perform_timing_ ) { gt_timer_.stop(); }

        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(complex_image, debug_folder_ + "raw_data_recon_Result_complex_image");
        if ( !debug_folder_.empty() ) gt_exporter_.export_array_complex(full_kspace, debug_folder_ + "raw_data_recon_Result_full_kspace");
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
        size_t e1, n, s;
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

        if ( !debug_folder_.empty() ) gt_exporter_.export_array(binning_obj_.time_stamp_, debug_folder_+"binning_obj_time_stamp");

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

            if ( !debug_folder_.empty() ) gt_exporter_.export_array(time_stamp, debug_folder_+"time_stamp" + os.str());
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_stamp, debug_folder_+"cpt_time_stamp" + os.str());

            float meanRR(0);

            if ( !hasDescending )
            {
                if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_stamp, debug_folder_+"cpt_time_stamp" + os.str());

                this->process_time_stamps(time_stamp, cpt_time_stamp, 
                                cpt_time_ratio, phs_time_stamp, phs_cpt_time_stamp, phs_cpt_time_ratio, 
                                ind_heart_beat, startingHB, endingHB, meanRR);
            }
            else
            {
                GDEBUG_STREAM("CmrKSpaceBinning -- Alternating sampling strategy detected ... ");

                if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_stamp_flipped, debug_folder_+"cpt_time_stamp_flipped" + os.str());

                hoNDArray<float> cpt_time_ratio_tmp(E1, N);

                this->process_time_stamps(time_stamp, cpt_time_stamp_flipped, 
                                cpt_time_ratio_tmp, phs_time_stamp, phs_cpt_time_stamp, phs_cpt_time_ratio, 
                                ind_heart_beat, startingHB, endingHB, meanRR);

                if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_ratio_tmp, debug_folder_+"cpt_time_ratio_tmp" + os.str());

                cpt_time_stamp = cpt_time_stamp_flipped;
                cpt_time_ratio = cpt_time_ratio_tmp;

                for ( n=0; n<N; n++ )
                {
                    if ( !ascending[n] )
                    {
                        for ( e1=0; e1<E1; e1++ )
                        {
                            cpt_time_ratio(e1, n) = cpt_time_ratio_tmp(E1-1-e1, n);
                        }
                    }
                }
            }

            if ( !debug_folder_.empty() ) gt_exporter_.export_array(cpt_time_ratio, debug_folder_+"cpt_time_ratio" + os.str());
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(phs_time_stamp, debug_folder_+"phs_time_stamp" + os.str());
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(phs_cpt_time_stamp, debug_folder_+"phs_cpt_time_stamp" + os.str());
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(phs_cpt_time_ratio, debug_folder_+"phs_cpt_time_ratio" + os.str());

            binning_obj_.starting_heart_beat_[s] = startingHB;
            binning_obj_.ending_heart_beat_[s] = endingHB;
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

        Gadgetron::hoImageRegContainer2DRegistration<T, float, 2, 2> reg;

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

            if ( !debug_folder_.empty() ) gt_exporter_.export_array(mag, debug_folder_+"navigator_mag" + os.str());

            // detect key frame
            size_t key_frame;
            Gadgetron::find_key_frame_SSD_2DT(mag, key_frame);
            GDEBUG_STREAM("Find key frame " << key_frame << " for S " << s);

            // perform motion correction
            Gadgetron:: perform_moco_fixed_key_frame_2DT(mag, key_frame, respiratory_navigator_moco_reg_strength_, respiratory_navigator_moco_iters_, false, reg);

            // get the moco results
            hoNDArray<T> magMoCo(mag);
            hoNDArray<float> dx, dy;

            reg.warped_container_.to_NDArray(0, magMoCo);
            reg.deformation_field_[0].to_NDArray(0, dx);
            reg.deformation_field_[1].to_NDArray(0, dy);

            if ( !debug_folder_.empty() ) gt_exporter_.export_array(magMoCo, debug_folder_+"navigator_mag_moco" + os.str());
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(dx, debug_folder_+"navigator_mag_moco_dx" + os.str());
            if ( !debug_folder_.empty() ) gt_exporter_.export_array(dy, debug_folder_+"navigator_mag_moco_dy" + os.str());

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
            hoNDArray<float> deform_patch(patch_size_RO, patch_size_E1);

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
                            memcpy(deform_patch.begin()+pl*patch_size_RO, dx.begin()+offset, sizeof(float)*patch_size_RO);
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
                            memcpy(deform_patch.begin()+pl*patch_size_RO, dy.begin()+offset, sizeof(float)*patch_size_RO);
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

            if ( !debug_folder_.empty() ) gt_exporter_.export_array(navigator_s, debug_folder_ + ""  + os.str());
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

                float mean_resp, var_resp;
                this->compute_metrics_navigator_heart_beat(s, HB, mean_resp, var_resp);

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

        long long n, lin;

        for ( n=0; n<N; n++ )
        {
            float v_first = -1;
            for ( lin=0; lin<E1; lin++ )
            {
                if ( time_stamp(lin, n) > -1 )
                {
                    v_first = time_stamp(lin, n);
                    break;
                }
            }

            float v_last = -1;
            for ( lin=E1-1; lin>=0; lin-- )
            {
                if ( time_stamp(lin, n) > -1 )
                {
                    v_last = time_stamp(lin, n);
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
                for ( lin=E1-1; lin>=0; lin-- )
                {
                    if ( cpt_time_stamp(lin, n) > -1 )
                    {
                        cpt_time_stamp_flipped(E1-1-lin, n) = cpt_time_stamp(lin, n);
                    }
                }
            }
            else
            {
                for ( lin=0; lin<E1; lin++ )
                {
                    if ( cpt_time_stamp(lin, n) > -1 )
                    {
                        cpt_time_stamp_flipped(lin, n) = cpt_time_stamp(lin, n);
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

        Gadgetron::correct_time_stamp_with_fitting(time_stamp);

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
        Gadgetron::correct_heart_beat_time_stamp_with_fitting(cpt_time_stamp, indHeartBeat, start_e1_hb, end_e1_hb, start_n_hb, end_n_hb);

        // --------------------------------------------------------
        // fill per phase time stamp  and phase cpt time stamp
        // --------------------------------------------------------
        Gadgetron::compute_phase_time_stamp(time_stamp, cpt_time_stamp, phs_time_stamp, phs_cpt_time_stamp);

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

                cpt_time_ratio(e1, n) = (cpt_time_stamp(e1, n) - minCPT) / (maxCPT-minCPT);
            }
        }

        // --------------------------------------------------------
        // fill the phase cpt time ratio
        // --------------------------------------------------------
        for ( n=0; n<N; n++ )
        {
            std::vector<float> buf(E1);
            for ( e1=0; e1<E1; e1++ )
            {
                buf[e1] = cpt_time_ratio(e1, n);
            }

            std::sort(buf.begin(), buf.end());
            phs_cpt_time_ratio(n, 0) = (float)(0.5*(buf[E1/2-1]+buf[E1/2]));
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
            if ( v < min_time_stamp ) min_time_stamp = v;
            if ( v > max_time_stamp ) max_time_stamp = v;
        }

        RRInterval = max_time_stamp - min_time_stamp;
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::compute_RRInterval() ... ");
    }
}

template <typename T> 
void CmrKSpaceBinning<T>::compute_metrics_navigator_heart_beat(size_t s, size_t HB, float& mean_resp, float& var_resp)
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
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in CmrKSpaceBinning<T>::compute_metrics_navigator_heart_beat() ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCMR CmrKSpaceBinning< float >;

}
