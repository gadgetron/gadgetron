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
CmrKSpaceBinning<T>::CmrKSpaceBinning()
{
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

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCMR CmrKSpaceBinning< float >;
template class EXPORTCMR CmrKSpaceBinning< double >;

}
