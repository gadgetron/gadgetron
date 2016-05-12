/** \file   cmr_kspace_binning.cpp
    \brief  Implement kspace binning recon for 2D acquisition
    \author Hui Xue
*/

#include "cmr_kspace_binning.h"
#include "hoNDFFT.h"
#include "mri_core_grappa.h"
#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_coil_map_estimation.h"

#include "cmr_time_stamp.h"

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
            startingHB[ind].second = start_n_hb[0];

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

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCMR CmrKSpaceBinning< float >;
template class EXPORTCMR CmrKSpaceBinning< double >;

}
