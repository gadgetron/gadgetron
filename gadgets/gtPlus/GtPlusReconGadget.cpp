
#include "GtPlusReconGadget.h"
#include "GtPlusGadgetOpenMP.h"
#include "gadgetron_paths.h"
#include <iomanip>
#include "CloudBus.h"

#include "mri_core_kspace_filter.h"

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

    GtPlusReconGadget::GtPlusReconGadget()
    {
        image_series_ = 100;

        min_intensity_value_ = 64;
        max_intensity_value_ = 4095;

        max_intensity_value_US_ = 2048;

        scalingFactor_ = -1;
        scalingFactor_gfactor_ = 100;
        scalingFactor_wrap_around_map_ = 1000;
        scalingFactor_snr_image_ = 10;
        scalingFactor_std_map_ = 1000;

        start_frame_for_std_map_ = 5;

        use_constant_scalingFactor_ = false;

        timeStampResolution_ = 0.0025f;

        aSpacing_[0] = 2.0;
        aSpacing_[1] = 2.0;
        aSpacing_[2] = 6.0;
        aSpacing_[3] = 1.0;
        aSpacing_[4] = 1.0;
        aSpacing_[5] = 1.0;

        reconE1_ = 1;
        reconE2_ = 1;

        processed_called_times_ = 0;

        thread_number_ratio_ = 0;

        kSpaceMaxAcqE2No_ = 0;

        min_E1_ = 0;
        max_E1_ = 0;
        center_E1_ = 0;

        min_E2_ = 0;
        max_E2_ = 0;
        center_E2_ = 0;

        filterRO_type_ = ISMRMRD_FILTER_GAUSSIAN;
        filterRO_sigma_ = 1.5;
        filterRO_width_ = 0.15;

        filterE1_type_ = ISMRMRD_FILTER_GAUSSIAN;
        filterE1_sigma_ = 1.5;
        filterE1_width_ = 0.15;

        filterE2_type_ = ISMRMRD_FILTER_GAUSSIAN;
        filterE2_sigma_ = 1.5;
        filterE2_width_ = 0.15;

        filterRO_ref_type_ = ISMRMRD_FILTER_HANNING;
        filterRO_ref_sigma_ = 1.5;
        filterRO_ref_width_ = 0.15;

        filterE1_ref_type_ = ISMRMRD_FILTER_HANNING;
        filterE1_ref_sigma_ = 1.5;
        filterE1_ref_width_ = 0.15;

        filterE2_ref_type_ = ISMRMRD_FILTER_HANNING;
        filterE2_ref_sigma_ = 1.5;
        filterE2_ref_width_ = 0.15;

        filterRO_pf_type_ = ISMRMRD_FILTER_HANNING;
        filterRO_pf_sigma_ = 1.5;
        filterRO_pf_width_ = 0.15;
        filterRO_pf_densityComp_ = false;

        filterE1_pf_type_ = ISMRMRD_FILTER_HANNING;
        filterE1_pf_sigma_ = 1.5;
        filterE1_pf_width_ = 0.15;
        filterE1_pf_densityComp_ = false;

        filterE2_pf_type_ = ISMRMRD_FILTER_HANNING;
        filterE2_pf_sigma_ = 1.5;
        filterE2_pf_width_ = 0.15;
        filterE2_pf_densityComp_ = false;

        recon_res_second_required_ = false;

        send_out_recon_ = true;
        send_out_recon_second_ = true;

        debugFolder_ = "DebugOutput";
        debugFolder2_ = debugFolder_;

        performTiming_ = true;

        verboseMode_ = false;

        CloudComputing_ = false;
        CloudSize_ = 0;

        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);

        Gadgetron::prepOpenMP();
    }

    GtPlusReconGadget::~GtPlusReconGadget()
    {
    }

    bool GtPlusReconGadget::readParameters()
    {
        try
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "------> GtPlusReconGadget parameters <------");

            min_intensity_value_ = min_intensity_value.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "min_intensity_value_ is " << min_intensity_value_);

            max_intensity_value_ = max_intensity_value.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "max_intensity_value_ is " << max_intensity_value_);

            scalingFactor_ = scalingFactor.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_ is " << scalingFactor_);

            scalingFactor_gfactor_ = scalingFactor_gfactor.value();
            if ( scalingFactor_gfactor_ == 0 ) scalingFactor_gfactor_ = 100;
            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_gfactor_ is " << scalingFactor_gfactor_);

            scalingFactor_wrap_around_map_ = scalingFactor_wrap_around_map.value();
            if ( scalingFactor_wrap_around_map_ == 0 ) scalingFactor_wrap_around_map_ = 1000;
            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_wrap_around_map_ is " << scalingFactor_wrap_around_map_);

            scalingFactor_snr_image_ = scalingFactor_snr_image.value();
            if ( scalingFactor_snr_image_ == 0 ) scalingFactor_snr_image_ = 10;
            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_snr_image_ is " << scalingFactor_snr_image_);

            scalingFactor_std_map_ = scalingFactor_std_map.value();
            if ( scalingFactor_std_map_ == 0 ) scalingFactor_std_map_ = 1000;
            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_std_map_ is " << scalingFactor_std_map_);

            start_frame_for_std_map_ = start_frame_for_std_map.value();
            if ( start_frame_for_std_map_ == 0 ) start_frame_for_std_map_ = 5;
            GDEBUG_CONDITION_STREAM(verboseMode_, "start_frame_for_std_map_ is " << start_frame_for_std_map_);

            use_constant_scalingFactor_ = use_constant_scalingFactor.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "use_constant_scalingFactor_ is " << use_constant_scalingFactor_);

            debugFolder_ = debugFolder.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "debugFolder_ is " << debugFolder_);

            debugFolder2_ = debugFolder2.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "debugFolder2_ is " << debugFolder2_);

            timeStampResolution_ = timeStampResolution.value();
            if ( timeStampResolution_ < FLT_EPSILON ) timeStampResolution_ = 0.0025f;
            GDEBUG_CONDITION_STREAM(verboseMode_, "timeStampResolution_ is " << timeStampResolution_);

            send_out_recon_ = send_out_recon.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "send_out_recon_ is " << send_out_recon_);

            send_out_recon_second_ = send_out_recon_second.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "send_out_recon_second_ is " << send_out_recon_second_);

            performTiming_ = performTiming.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "performTiming_ is " << performTiming_);

            // kspace filter parameters
            filterRO_type_ = Gadgetron::get_kspace_filter_type(filterRO.value());
            filterRO_sigma_ = filterRO_sigma.value();
            filterRO_width_ = filterRO_width.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_type_ is " << filterRO.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_sigma_ is " << filterRO_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_width_ is " << filterRO_width_);

            filterE1_type_ = Gadgetron::get_kspace_filter_type(filterE1.value());
            filterE1_sigma_ = filterE1_sigma.value();
            filterE1_width_ = filterE1_width.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_type_ is " << filterE1.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_sigma_ is " << filterE1_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_width_ is " << filterE1_width_);

            filterE2_type_ = Gadgetron::get_kspace_filter_type(filterE2.value());
            filterE2_sigma_ = filterE2_sigma.value();
            filterE2_width_ = filterE2_width.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_type_ is " << filterE2.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_sigma_ is " << filterE2_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_width_ is " << filterE2_width_);

            filterRO_ref_type_ = Gadgetron::get_kspace_filter_type(filterRefRO.value());
            filterRO_ref_sigma_ = filterRefRO_sigma.value();
            filterRO_ref_width_ = filterRefRO_width.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_ref_type_ is " << filterRefRO.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_ref_sigma_ is " << filterRO_ref_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_ref_width_ is " << filterRO_ref_width_);

            filterE1_ref_type_ = Gadgetron::get_kspace_filter_type(filterRefE1.value());
            filterE1_ref_sigma_ = filterRefE1_sigma.value();
            filterE1_ref_width_ = filterRefE1_width.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_ref_type_ is " << filterRefE1.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_ref_sigma_ is " << filterE1_ref_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_ref_width_ is " << filterE1_ref_width_);

            filterE2_ref_type_ = Gadgetron::get_kspace_filter_type(filterRefE2.value());
            filterE2_ref_sigma_ = filterRefE2_sigma.value();
            filterE2_ref_width_ = filterRefE2_width.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_ref_type_ is " << filterRefE2.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_ref_sigma_ is " << filterE2_ref_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_ref_width_ is " << filterE2_ref_width_);

            filterRO_pf_type_ = Gadgetron::get_kspace_filter_type(filterPartialFourierRO.value());
            filterRO_pf_sigma_ = filterPartialFourierRO_sigma.value();
            filterRO_pf_width_ = filterPartialFourierRO_width.value();
            filterRO_pf_densityComp_ = filterPartialFourierRO_densityComp.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_type_ is " << filterPartialFourierRO.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_sigma_ is " << filterRO_pf_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_width_ is " << filterRO_pf_width_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_densityComp_ is " << filterRO_pf_densityComp_);

            filterE1_pf_type_ = Gadgetron::get_kspace_filter_type(filterPartialFourierE1.value());
            filterE1_pf_sigma_ = filterPartialFourierE1_sigma.value();
            filterE1_pf_width_ = filterPartialFourierE1_width.value();
            filterE1_pf_densityComp_ = filterPartialFourierE1_densityComp.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_type_ is " << filterPartialFourierE1.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_sigma_ is " << filterE1_pf_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_width_ is " << filterE1_pf_width_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_densityComp_ is " << filterE1_pf_densityComp_);

            filterE2_pf_type_ = Gadgetron::get_kspace_filter_type(filterPartialFourierE2.value());
            filterE2_pf_sigma_ = filterPartialFourierE2_sigma.value();
            filterE2_pf_width_ = filterPartialFourierE2_width.value();
            filterE2_pf_densityComp_ = filterPartialFourierE2_densityComp.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_type_ is " << filterPartialFourierE2.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_sigma_ is " << filterE2_pf_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_width_ is " << filterE2_pf_width_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_densityComp_ is " << filterE2_pf_densityComp_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            CloudComputing_ = CloudComputing.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "CloudComputing_ is " << CloudComputing_);

            cloud_node_file_ = cloudNodeFile.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "cloud_node_file_ is " << cloud_node_file_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "CloudNodeXMLConfiguration is " << CloudNodeXMLConfiguration.value());

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            thread_number_ratio_ = 1;
            if ( thread_number_ratio_>1 || thread_number_ratio_<0 ) thread_number_ratio_ = 0;
            GDEBUG_CONDITION_STREAM(verboseMode_, "thread_number_ratio_ is " << thread_number_ratio_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "==================================================================");

            GDEBUG_CONDITION_STREAM(verboseMode_, "------> GtPlus recon parameters <------");

            workOrderPara_.upstream_coil_compression_ = upstream_coil_compression.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "upstream_coil_compression_ is " << workOrderPara_.upstream_coil_compression_);

            workOrderPara_.upstream_coil_compression_thres_ = upstream_coil_compression_thres.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "upstream_coil_compression_thres_ is " << workOrderPara_.upstream_coil_compression_thres_);

            workOrderPara_.upstream_coil_compression_num_modesKept_ = upstream_coil_compression_num_modesKept.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "upstream_coil_compression_num_modesKept_ is " << workOrderPara_.upstream_coil_compression_num_modesKept_);

            workOrderPara_.downstream_coil_compression_ = downstream_coil_compression.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "downstream_coil_compression_ is " << workOrderPara_.downstream_coil_compression_);

            workOrderPara_.coil_compression_thres_ = coil_compression_thres.value();

            if (workOrderPara_.upstream_coil_compression_ && (workOrderPara_.upstream_coil_compression_thres_>0) && (workOrderPara_.coil_compression_thres_ > workOrderPara_.upstream_coil_compression_thres_))
                workOrderPara_.coil_compression_thres_ = workOrderPara_.upstream_coil_compression_thres_;

            GDEBUG_CONDITION_STREAM(verboseMode_, "coil_compression_thres_ is " << workOrderPara_.coil_compression_thres_);

            workOrderPara_.coil_compression_num_modesKept_ = coil_compression_num_modesKept.value();

            if (workOrderPara_.upstream_coil_compression_ && (workOrderPara_.upstream_coil_compression_num_modesKept_>0) && (workOrderPara_.coil_compression_num_modesKept_ > workOrderPara_.upstream_coil_compression_num_modesKept_))
                workOrderPara_.coil_compression_num_modesKept_ = workOrderPara_.upstream_coil_compression_num_modesKept_;

            GDEBUG_CONDITION_STREAM(verboseMode_, "coil_compression_num_modesKept_ is " << workOrderPara_.coil_compression_num_modesKept_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.coil_map_algorithm_ = gtPlus_util_.getISMRMRDCoilMapAlgoFromName(coil_map_algorithm.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "coil_map_algorithm_ is " << coil_map_algorithm.value());

            workOrderPara_.csm_kSize_ = csm_kSize.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_kSize_ is " << workOrderPara_.csm_kSize_);

            workOrderPara_.csm_powermethod_num_ = csm_powermethod_num.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_powermethod_num_ is " << workOrderPara_.csm_powermethod_num_);

            workOrderPara_.csm_true_3D_ = csm_true_3D.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_true_3D_ is " << workOrderPara_.csm_true_3D_);

            workOrderPara_.csm_iter_num_ = csm_iter_num.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_iter_num_ is " << workOrderPara_.csm_iter_num_);

            workOrderPara_.csm_iter_thres_ = csm_iter_thres.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_iter_thres_ is " << workOrderPara_.csm_iter_thres_);

            workOrderPara_.csm_use_gpu_ = false;
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_use_gpu_ is " << workOrderPara_.csm_use_gpu_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.recon_algorithm_ = gtPlus_util_.getISMRMRDReconAlgoFromName(recon_algorithm.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "recon_algorithm_ is " << recon_algorithm.value());

            workOrderPara_.recon_auto_parameters_ = recon_auto_parameters.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "recon_auto_parameters_ is " << workOrderPara_.recon_auto_parameters_);

            workOrderPara_.gfactor_needed_ = gfactor_needed.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "gfactor_needed_ is " << workOrderPara_.gfactor_needed_);

            workOrderPara_.wrap_around_map_needed_ = wrap_around_map_needed.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "wrap_around_map_needed_ is " << workOrderPara_.wrap_around_map_needed_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.grappa_kSize_RO_ = grappa_kSize_RO.value();
            workOrderPara_.grappa_kSize_E1_ = grappa_kSize_E1.value();
            workOrderPara_.grappa_kSize_E2_ = grappa_kSize_E2.value();
            workOrderPara_.grappa_reg_lamda_ = grappa_reg_lamda.value();
            workOrderPara_.grappa_calib_over_determine_ratio_ = grappa_calib_over_determine_ratio.value();
            workOrderPara_.grappa_use_gpu_ = false;

            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_kSize_RO_ is " << workOrderPara_.grappa_kSize_RO_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_kSize_E1_ is " << workOrderPara_.grappa_kSize_E1_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_kSize_E2_ is " << workOrderPara_.grappa_kSize_E2_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_reg_lamda_ is " << workOrderPara_.grappa_reg_lamda_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_calib_over_determine_ratio_ is " << workOrderPara_.grappa_calib_over_determine_ratio_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_use_gpu_ is " << workOrderPara_.grappa_use_gpu_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.spirit_kSize_RO_ = spirit_kSize_RO.value();
            if ( workOrderPara_.spirit_kSize_RO_ == 0 ) workOrderPara_.spirit_kSize_RO_ = 7;

            workOrderPara_.spirit_kSize_E1_ = spirit_kSize_E1.value();
            if ( workOrderPara_.spirit_kSize_E1_ == 0 ) workOrderPara_.spirit_kSize_E1_ = 7;

            workOrderPara_.spirit_kSize_E2_ = spirit_kSize_E2.value();
            if ( workOrderPara_.spirit_kSize_E2_ == 0 ) workOrderPara_.spirit_kSize_E2_ = 5;

            workOrderPara_.spirit_reg_lamda_ = spirit_reg_lamda.value();
            if ( workOrderPara_.spirit_reg_lamda_ < FLT_EPSILON ) workOrderPara_.spirit_reg_lamda_ = 0.005;

            workOrderPara_.spirit_use_gpu_ = false;
            workOrderPara_.spirit_calib_over_determine_ratio_ = spirit_calib_over_determine_ratio.value();
            workOrderPara_.spirit_solve_symmetric_ = spirit_solve_symmetric.value();

            workOrderPara_.spirit_iter_max_ = spirit_iter_max.value();
            if ( workOrderPara_.spirit_iter_max_ == 0 ) workOrderPara_.spirit_iter_max_ = 100;

            workOrderPara_.spirit_iter_thres_ = spirit_iter_thres.value();
            if ( workOrderPara_.spirit_iter_thres_ < FLT_EPSILON ) workOrderPara_.spirit_iter_thres_ = 0.0015;

            workOrderPara_.spirit_print_iter_ = spirit_print_iter.value();

            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_kSize_RO_ is " << workOrderPara_.spirit_kSize_RO_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_kSize_E1_ is " << workOrderPara_.spirit_kSize_E1_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_kSize_E2_ is " << workOrderPara_.spirit_kSize_E2_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_reg_lamda_ is " << workOrderPara_.spirit_reg_lamda_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_use_gpu_ is " << workOrderPara_.spirit_use_gpu_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_calib_over_determine_ratio_ is " << workOrderPara_.spirit_calib_over_determine_ratio_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_solve_symmetric_ is " << workOrderPara_.spirit_solve_symmetric_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_iter_max_ is " << workOrderPara_.spirit_iter_max_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_iter_thres_ is " << workOrderPara_.spirit_iter_thres_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_print_iter_ is " << workOrderPara_.spirit_print_iter_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.spirit_perform_linear_ = spirit_perform_linear.value();
            workOrderPara_.spirit_perform_nonlinear_ = spirit_perform_nonlinear.value();
            workOrderPara_.spirit_parallel_imaging_lamda_ = spirit_parallel_imaging_lamda.value();
            workOrderPara_.spirit_image_reg_lamda_ = spirit_image_reg_lamda.value();
            workOrderPara_.spirit_data_fidelity_lamda_ = spirit_data_fidelity_lamda.value();
            workOrderPara_.spirit_ncg_iter_max_ = spirit_ncg_iter_max.value();
            workOrderPara_.spirit_ncg_iter_thres_ = spirit_ncg_iter_thres.value();
            workOrderPara_.spirit_ncg_print_iter_ = spirit_ncg_print_iter.value();
            // spirit_ncg_scale_factor_ is computed from the data

            workOrderPara_.spirit_use_coil_sen_map_ = spirit_use_coil_sen_map.value();
            workOrderPara_.spirit_use_moco_enhancement_ = spirit_use_moco_enhancement.value();
            workOrderPara_.spirit_recon_moco_images_ = spirit_recon_moco_images.value();
            workOrderPara_.spirit_RO_enhancement_ratio_ = spirit_RO_enhancement_ratio.value();
            workOrderPara_.spirit_E1_enhancement_ratio_ = spirit_E1_enhancement_ratio.value();
            workOrderPara_.spirit_E2_enhancement_ratio_ = spirit_E2_enhancement_ratio.value();
            workOrderPara_.spirit_temporal_enhancement_ratio_ = spirit_temporal_enhancement_ratio.value();
            workOrderPara_.spirit_2D_scale_per_chunk_ = spirit_2D_scale_per_chunk.value();
            workOrderPara_.spirit_3D_scale_per_chunk_ = spirit_3D_scale_per_chunk.value();

            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_perform_linear_ is " << workOrderPara_.spirit_perform_linear_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_perform_nonlinear_ is " << workOrderPara_.spirit_perform_nonlinear_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_parallel_imaging_lamda_ is " << workOrderPara_.spirit_parallel_imaging_lamda_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_image_reg_lamda_ is " << workOrderPara_.spirit_image_reg_lamda_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_data_fidelity_lamda_ is " << workOrderPara_.spirit_data_fidelity_lamda_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_ncg_iter_max_ is " << workOrderPara_.spirit_ncg_iter_max_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_ncg_iter_thres_ is " << workOrderPara_.spirit_ncg_iter_thres_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_ncg_print_iter_ is " << workOrderPara_.spirit_ncg_print_iter_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_use_coil_sen_map_ is " << workOrderPara_.spirit_use_coil_sen_map_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_use_moco_enhancement_ is " << workOrderPara_.spirit_use_moco_enhancement_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_recon_moco_images_ is " << workOrderPara_.spirit_recon_moco_images_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_RO_enhancement_ratio_ is " << workOrderPara_.spirit_RO_enhancement_ratio_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_E1_enhancement_ratio_ is " << workOrderPara_.spirit_E1_enhancement_ratio_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_E2_enhancement_ratio_ is " << workOrderPara_.spirit_E2_enhancement_ratio_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_temporal_enhancement_ratio_ is " << workOrderPara_.spirit_temporal_enhancement_ratio_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_2D_scale_per_chunk_ is " << workOrderPara_.spirit_2D_scale_per_chunk_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "spirit_3D_scale_per_chunk_ is " << workOrderPara_.spirit_3D_scale_per_chunk_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.retro_gated_interp_method_ = gtPlus_util_.getISMRMRDRetroGatingInterpFromName(retro_gated_interp_method.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "retro_gated_interp_method_ is " << retro_gated_interp_method.value());

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.job_split_by_S_ = job_split_by_S.value();
            workOrderPara_.job_num_of_N_ = job_num_of_N.value();
            workOrderPara_.job_max_Megabytes_ = job_max_Megabytes.value();
            workOrderPara_.job_overlap_ = job_overlap.value();
            workOrderPara_.job_perform_on_control_node_ = job_perform_on_control_node.value();

            GDEBUG_CONDITION_STREAM(verboseMode_, "job_split_by_S_ is " << workOrderPara_.job_split_by_S_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "job_num_of_N_ is " << workOrderPara_.job_num_of_N_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "job_max_Megabytes_ is " << workOrderPara_.job_max_Megabytes_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "job_overlap_ is " << workOrderPara_.job_overlap_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "job_perform_on_control_node_ is " << workOrderPara_.job_perform_on_control_node_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.partialFourier_algo_ = gtPlus_util_.getISMRMRDPartialFourierReconAlgoFromName(partialFourier_algo.value());
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_algo_ is " << partialFourier_algo.value());

            workOrderPara_.partialFourier_homodyne_iters_ = partialFourier_homodyne_iters.value();
            workOrderPara_.partialFourier_homodyne_thres_ = partialFourier_homodyne_thres.value();
            workOrderPara_.partialFourier_homodyne_densityComp_ = partialFourier_homodyne_densityComp.value();

            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_homodyne_iters_ is " << workOrderPara_.partialFourier_homodyne_iters_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_homodyne_thres_ is " << workOrderPara_.partialFourier_homodyne_thres_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_homodyne_densityComp_ is " << workOrderPara_.partialFourier_homodyne_densityComp_);

            workOrderPara_.partialFourier_POCS_iters_ = partialFourier_POCS_iters.value();
            workOrderPara_.partialFourier_POCS_thres_ = partialFourier_POCS_thres.value();
            workOrderPara_.partialFourier_POCS_transitBand_ = partialFourier_POCS_transitBand.value();
            workOrderPara_.partialFourier_POCS_transitBand_E2_ = partialFourier_POCS_transitBand_E2.value();

            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_POCS_iters_ is " << workOrderPara_.partialFourier_POCS_iters_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_POCS_thres_ is " << workOrderPara_.partialFourier_POCS_thres_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_POCS_transitBand_ is " << workOrderPara_.partialFourier_POCS_transitBand_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_POCS_transitBand_ is " << workOrderPara_.partialFourier_POCS_transitBand_E2_);

            workOrderPara_.partialFourier_FengHuang_kSize_RO_ = partialFourier_FengHuang_kSize_RO.value();
            workOrderPara_.partialFourier_FengHuang_kSize_E1_ = partialFourier_FengHuang_kSize_E1.value();
            workOrderPara_.partialFourier_FengHuang_kSize_E2_ = partialFourier_FengHuang_kSize_E2.value();
            workOrderPara_.partialFourier_FengHuang_thresReg_ = partialFourier_FengHuang_thresReg.value();
            workOrderPara_.partialFourier_FengHuang_sameKernel_allN_ = partialFourier_FengHuang_sameKernel_allN.value();
            workOrderPara_.partialFourier_FengHuang_transitBand_ = partialFourier_FengHuang_transitBand.value();
            workOrderPara_.partialFourier_FengHuang_transitBand_E2_ = partialFourier_FengHuang_transitBand_E2.value();

            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_kSize_RO_ is " << workOrderPara_.partialFourier_FengHuang_kSize_RO_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_kSize_E1_ is " << workOrderPara_.partialFourier_FengHuang_kSize_E1_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_kSize_E2_ is " << workOrderPara_.partialFourier_FengHuang_kSize_E2_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_thresReg_ is " << workOrderPara_.partialFourier_FengHuang_thresReg_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_sameKernel_allN_ is " << workOrderPara_.partialFourier_FengHuang_sameKernel_allN_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_transitBand_ is " << workOrderPara_.partialFourier_FengHuang_transitBand_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_transitBand_E2_ is " << workOrderPara_.partialFourier_FengHuang_transitBand_E2_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            recon_kspace_needed_ = recon_kspace_needed.value();
            GDEBUG_CONDITION_STREAM(verboseMode_, "recon_kspace_needed_ is " << recon_kspace_needed_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");
        }
        catch(...)
        {
            GERROR_STREAM("Errors in GtPlusReconGadget::readParameters() ... ");
            return false;
        }

        return true;
    }

    bool GtPlusReconGadget::parseGTCloudNodeFile(const std::string& filename, CloudType& gtCloud)
    {

        bool has_cloud_node_xml_configuration = true;

        if (this->using_cloudbus.value() && has_cloud_node_xml_configuration) {
            std::vector<GadgetronNodeInfo> nodes;
            CloudBus::instance()->get_node_info(nodes);

            if (nodes.size()>0)
            {
                gtCloud.resize(nodes.size());

                unsigned int n;
                for ( n=0; n<nodes.size(); n++ )
                {
                    std::stringstream ss;
                    gtCloud[n].get<0>() = nodes[n].address;
                    ss << nodes[n].port;
                    gtCloud[n].get<1>() = ss.str();
                    gtCloud[n].get<2>() = CloudNodeXMLConfiguration.value();
                    gtCloud[n].get<3>() = nodes[n].compute_capability;

                    GDEBUG_CONDITION_STREAM(verboseMode_, "Gadget Node " << n << " : " << gt_cloud_[n]);
                }

                return true; //We will leave the function here
            }
            else
            {
                GDEBUG_STREAM("Cloud bus cannot find any nodes; using the node text file instead ... ");
            }
        }

        std::string nodeFileName = get_gadgetron_home();
        nodeFileName.append("/share/gadgetron/config/gtCloud/");
        nodeFileName.append(filename);
        GDEBUG_CONDITION_STREAM(verboseMode_, "Cloud node file name is " << nodeFileName);

        std::ifstream fs(nodeFileName.c_str(), std::ios::in);
        if (!fs.is_open()) 
        {
            GWARN_STREAM("Cannot open GT CloudNodeFile; use the local setting instead ... ");
            return false;
        }

        // control node hostname
        std::string controlNode;
        fs >> controlNode;

        std::string portControlNode;
        fs >> portControlNode;

        // number of GadgetLevel nodes
        unsigned int num;
        fs >> num;

        gtCloud.resize(num);

        unsigned int n;
        for ( n=0; n<num; n++ )
        {
            std::string gadgetNode;
            fs >> gadgetNode;

            std::string portGadgetNode;
            fs >> portGadgetNode;

            std::string xmlGadgetNode;
            fs >> xmlGadgetNode;

            unsigned int computingPowerIndex;
            fs >> computingPowerIndex;

            gtCloud[n].get<0>() = gadgetNode;
            gtCloud[n].get<1>() = portGadgetNode;
            gtCloud[n].get<2>() = xmlGadgetNode;
            gtCloud[n].get<3>() = computingPowerIndex;

            GDEBUG_CONDITION_STREAM(verboseMode_, "Gadget Node " << n << " : " << gt_cloud_[n]);
        }

        fs.close();

        return true;
    }

    int GtPlusReconGadget::process_config(ACE_Message_Block* mb)
    {
        // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave]
        //   0  1  2   3    4   5    6     7  8   9   10

        verboseMode_ = verboseMode.value();

        // read parameters from xml
        image_series_ = image_series.value();

        // read in parameters from the xml
        GADGET_CHECK_RETURN(this->readParameters(), GADGET_FAIL);

        // check whether the second set of recon results is required
        recon_res_second_required_ = false;

        ISMRMRD::IsmrmrdHeader h;
        try {
            deserialize(mb->rd_ptr(),h);
        } catch (...) {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!h.acquisitionSystemInformation) {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }
        num_acq_channels_ = h.acquisitionSystemInformation->receiverChannels;

        GDEBUG_CONDITION_STREAM(verboseMode_, "Number of acquisition channels : " << num_acq_channels_);

        if (h.encoding.size() < 1 || h.encoding.size() > 2) {
            GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
            GDEBUG("This GtPlusReconGadget only supports one or two encoding spaces\n");
            return GADGET_FAIL;
        }

        // find out the encoding space 
        ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
        ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
        ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

        matrix_size_encoding_[0] = e_space.matrixSize.x;
        matrix_size_encoding_[1] = e_space.matrixSize.y;
        matrix_size_encoding_[2] = e_space.matrixSize.z;
        GDEBUG_CONDITION_STREAM(verboseMode_, "Encoding matrix size: " << matrix_size_encoding_[0] << " " << matrix_size_encoding_[1] << " " << matrix_size_encoding_[2]);

        field_of_view_encoding_[0] = e_space.fieldOfView_mm.x;
        field_of_view_encoding_[1] = e_space.fieldOfView_mm.y;
        field_of_view_encoding_[2] = e_space.fieldOfView_mm.z;
        GDEBUG_CONDITION_STREAM(verboseMode_, "Encoding field_of_view : " << field_of_view_encoding_[0] << " " << field_of_view_encoding_[1] << " " << field_of_view_encoding_[2]);

        // find the recon space
        matrix_size_recon_[0] = r_space.matrixSize.x;
        matrix_size_recon_[1] = r_space.matrixSize.y;
        matrix_size_recon_[2] = r_space.matrixSize.z;
        GDEBUG_CONDITION_STREAM(verboseMode_, "Recon matrix size : " << matrix_size_recon_[0] << " " << matrix_size_recon_[1] << " " << matrix_size_recon_[2]);

        field_of_view_recon_[0] = r_space.fieldOfView_mm.x;
        field_of_view_recon_[1] = r_space.fieldOfView_mm.y;
        field_of_view_recon_[2] = r_space.fieldOfView_mm.z;
        GDEBUG_CONDITION_STREAM(verboseMode_, "Recon field_of_view :  " << field_of_view_recon_[0] << " " << field_of_view_recon_[1] << " " << field_of_view_recon_[2]);

        if (e_limits.kspace_encoding_step_1.is_present())
        {
            min_E1_ = e_limits.kspace_encoding_step_1->minimum;
            max_E1_ = e_limits.kspace_encoding_step_1->maximum;
            center_E1_ = e_limits.kspace_encoding_step_1->center;

            GDEBUG_CONDITION_STREAM(verboseMode_, "Encoding limit E1 : " << min_E1_ << " " << max_E1_ << " " << center_E1_);
        }

        if (e_limits.kspace_encoding_step_2.is_present())
        {
            min_E2_ = e_limits.kspace_encoding_step_2->minimum;
            max_E2_ = e_limits.kspace_encoding_step_2->maximum;
            center_E2_ = e_limits.kspace_encoding_step_2->center;

            GDEBUG_CONDITION_STREAM(verboseMode_, "Encoding limit E2 : " << min_E2_ << " " << max_E2_ << " " << center_E2_);
        }

        // this gadget supports two encoding spaces only if the
        // second encoding space has the same field of view and resolution as the first
        // e.g. for FLASH PAT reference scans.
        if (h.encoding.size() == 2)
        {
            if (! ((h.encoding[0].reconSpace.matrixSize.x == h.encoding[1].reconSpace.matrixSize.x) && 
                (h.encoding[0].reconSpace.matrixSize.y == h.encoding[1].reconSpace.matrixSize.y) && 
                (h.encoding[0].reconSpace.matrixSize.z == h.encoding[1].reconSpace.matrixSize.z) && 
                (h.encoding[0].reconSpace.fieldOfView_mm.x == h.encoding[1].reconSpace.fieldOfView_mm.x) &&
                (h.encoding[0].reconSpace.fieldOfView_mm.y == h.encoding[1].reconSpace.fieldOfView_mm.y) &&
                (h.encoding[0].reconSpace.fieldOfView_mm.z == h.encoding[1].reconSpace.fieldOfView_mm.z)) )
            {
                GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
                GDEBUG("This GtPlusAccumulatorWorkOrderTriggerGadget only supports two encoding spaces with identical recon spaces.\n");
                return GADGET_FAIL;
            }
        }

        reconE1_ = matrix_size_recon_[1];
        GDEBUG_CONDITION_STREAM(verboseMode_, "reconE1_ is " << reconE1_);

        reconE2_ = matrix_size_recon_[2];
        GDEBUG_CONDITION_STREAM(verboseMode_, "reconE2_ is " << reconE2_);

        kSpaceMaxAcqE1No_ = matrix_size_encoding_[1]-1;
        GDEBUG_CONDITION_STREAM(verboseMode_, "kSpaceMaxAcqE1No_ is " << kSpaceMaxAcqE1No_);

        kSpaceMaxAcqE2No_ = matrix_size_encoding_[2]-1;
        GDEBUG_CONDITION_STREAM(verboseMode_, "kSpaceMaxAcqE2No_ is " << kSpaceMaxAcqE2No_);

        aSpacing_[0] = field_of_view_recon_[0]/matrix_size_recon_[0];
        aSpacing_[1] = field_of_view_recon_[1]/reconE1_;
        aSpacing_[2] = field_of_view_recon_[2]/reconE2_;

        gt_exporter_.setPixelSize(aSpacing_[0], aSpacing_[1], aSpacing_[2], aSpacing_[3], aSpacing_[4], aSpacing_[5]);

        //XUE-TODO: This is actually wrong. This assumes that you always zeropad, which is probably bad practice
        meas_max_idx_.kspace_encode_step_1 = (uint16_t)matrix_size_encoding_[1]-1;
        meas_max_idx_.set = (e_limits.set && (e_limits.set->maximum>0)) ? e_limits.set->maximum : 0;
        meas_max_idx_.phase = (e_limits.phase && (e_limits.phase->maximum>0)) ? e_limits.phase->maximum : 0;

        meas_max_idx_.kspace_encode_step_2 = (uint16_t)matrix_size_encoding_[2]-1; 

        meas_max_idx_.contrast = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum : 0;
        meas_max_idx_.slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;
        meas_max_idx_.repetition = e_limits.repetition ? e_limits.repetition->maximum : 0;
        meas_max_idx_.average = e_limits.average ? e_limits.average->maximum : 0;

        // combine all incoming segments
        meas_max_idx_.segment = 0;

        // find out the PAT mode
        if (!h.encoding[0].parallelImaging) {
            GDEBUG("Parallel Imaging section not found in header");
            return GADGET_FAIL;
        }

        ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;

        acceFactorE1_ = (long)(p_imaging.accelerationFactor.kspace_encoding_step_1);
        acceFactorE2_ = (long)(p_imaging.accelerationFactor.kspace_encoding_step_2);
        GDEBUG_CONDITION_STREAM(verboseMode_, "acceFactorE1 is " << acceFactorE1_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "acceFactorE2 is " << acceFactorE2_);

        std::string calib = *p_imaging.calibrationMode;

        bool separate = (calib.compare("separate") == 0);
        bool embedded = (calib.compare("embedded") == 0);
        bool external = (calib.compare("external") == 0);
        bool interleaved = (calib.compare("interleaved") == 0);
        bool other = (calib.compare("other") == 0);

        if ( separate )
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Colibration mode is separate");
        }
        else if ( embedded )
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Colibration mode is embedded");
        }
        else if ( interleaved )
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Colibration mode is interleaved");
        }
        else if ( external )
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Colibration mode is external");
        }
        else if ( other )
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Colibration mode is other");
        }

        //if ( other_ && acceFactorE1_==1 && acceFactorE2_==1 )
        //{
        //    GDEBUG_CONDITION_STREAM(verboseMode_, "Colibration mode is changed to ISMRMRD_interleaved");
        //    CalibMode_ = Gadgetron::ISMRMRD_interleaved;
        //    acceFactorE1_ = 2;
        //}

        CalibMode_ = Gadgetron::ISMRMRD_noacceleration;

        if ( interleaved )
        {
            CalibMode_ = Gadgetron::ISMRMRD_interleaved;

            if ( p_imaging.interleavingDimension )
            {
                if ( p_imaging.interleavingDimension->compare("phase") == 0 ) {
                    InterleaveDim_ = Gadgetron::DIM_Phase;
                } else if ( p_imaging.interleavingDimension->compare("repetition") == 0 ) {
                    InterleaveDim_ = Gadgetron::DIM_Repetition;
                } else if ( p_imaging.interleavingDimension->compare("average") == 0 ) {
                    InterleaveDim_ = Gadgetron::DIM_Average;
                } else if ( p_imaging.interleavingDimension->compare("contrast") == 0 ) {
                    InterleaveDim_ = Gadgetron::DIM_Contrast;
                } else if ( p_imaging.interleavingDimension->compare("other") == 0 ) {
                    InterleaveDim_ = Gadgetron::DIM_other1;
                } else {
                    GDEBUG("Unknown interleaving dimension. Bailing out");
                    return GADGET_FAIL;
                }
            }
        }
        else if ( embedded )
        {
            CalibMode_ = Gadgetron::ISMRMRD_embedded;
        }
        else if ( separate )
        {
            CalibMode_ = Gadgetron::ISMRMRD_separate;
        }
        else if ( external )
        {
            CalibMode_ = Gadgetron::ISMRMRD_external;
        }
        else if ( other )
        {
            CalibMode_ = Gadgetron::ISMRMRD_other;
        }

        // ---------------------------------------------------------------------------------------------------------
        // generate the destination folder
        if ( !debugFolder_.empty() )
        {
            Gadgetron::getDebugFolderPath(debugFolder_, debugFolder_fullPath_, verboseMode_);
        }
        else
        {
            GDEBUG_STREAM("GtPlusRecon, debugFolder is not set ...");
        }

        if ( !debugFolder2_.empty() )
        {
            Gadgetron::getDebugFolderPath(debugFolder2_, debugFolder2_fullPath_, verboseMode_);
        }
        else
        {
            GDEBUG_STREAM("GtPlusRecon, debugFolder2 is not set ...");
        }

        // ---------------------------------------------------------------------------------------------------------
        // set the maximal number of threads used
        if ( thread_number_ratio_>0 && thread_number_ratio_<1 )
        {
        }

        return GADGET_OK;
    }

    int GtPlusReconGadget::process(Gadgetron::GadgetContainerMessage< GtPlusGadgetImageArray >* m1, Gadgetron::GadgetContainerMessage< WorkOrderType > * m2)
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusReconGadget::process(...) starts ... ");

        processed_called_times_++;

        GtPlusGadgetImageArray* images = m1->getObjectPtr();

        boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->data_.get_dimensions();

        GDEBUG_CONDITION_STREAM(verboseMode_, "[Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave] = [" 
            << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " 
            << (*dims)[3] << " " << (*dims)[4] << " " << (*dims)[5] << " " 
            << (*dims)[6] << " " << (*dims)[7] << " " << (*dims)[8] << " " 
            << (*dims)[9] << " " << (*dims)[10] << "]");

        dimensions_ = *dims;

        GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusReconGadget::process(...) ends ... ");

        m1->release();
        return GADGET_OK;
    }

    size_t GtPlusReconGadget::computeSeriesImageNumber (ISMRMRD::ImageHeader& imheader, size_t nCHA, size_t cha, size_t nE2, size_t e2)
    {
        size_t nSET = meas_max_idx_.set+1;
        size_t nREP = meas_max_idx_.repetition+1;
        size_t nPHS = meas_max_idx_.phase+1;
        size_t nSLC = meas_max_idx_.slice+1;
        size_t nCON = meas_max_idx_.contrast+1;
        if ( nE2 == 0 ) nE2 = 1;

        size_t imageNum = imheader.average*nREP*nSET*nPHS*nCON*nSLC*nE2*nCHA 
            + imheader.repetition*nSET*nPHS*nCON*nSLC*nE2*nCHA 
            + imheader.set*nPHS*nCON*nSLC*nE2*nCHA 
            + imheader.phase*nCON*nSLC*nE2*nCHA 
            + imheader.contrast*nSLC*nE2*nCHA
            + imheader.slice*nE2*nCHA 
            + e2*nCHA 
            + cha 
            + 1;

        return imageNum;
    }

    bool GtPlusReconGadget::
        addPrePostZeros(int centreNo, int sampleNo, int& PrePostZeros)
    {
        // 1 : pre zeros
        // 2 : post zeros
        // 0 : no zeros
        PrePostZeros = 0;

        if ( sampleNo <= 1 )
            return true;

        if ( 2*centreNo == sampleNo )
        {
            PrePostZeros = 0;
        }

        if ( 2*centreNo < sampleNo )
        {
            PrePostZeros = 1;
        }

        if ( 2*centreNo > sampleNo )
        {
            PrePostZeros = 2;
        }

        return true;
    }

    bool GtPlusReconGadget::
        scalingImages(hoNDArray<ValueType>& res)
    {
        if ( scalingFactor_ < 0 && !use_constant_scalingFactor_ )
        {
            hoNDArray<float> mag(res.get_dimensions());
            Gadgetron::abs(res, mag);
            GADGET_CHECK_RETURN_FALSE(this->scalingMagnitude(mag));
        }

        scal((float)scalingFactor_, res);

        return true;
    }

    bool GtPlusReconGadget::
        scalingMagnitude(hoNDArray<float>& mag)
    {
        if ( scalingFactor_ < 0 && !use_constant_scalingFactor_ )
        {
            // perform the scaling to [0 max_inten_value_]
            size_t ind;
            float maxInten;

            size_t RO = mag.get_size(0);
            size_t E1 = mag.get_size(1);
            size_t num = mag.get_number_of_elements()/(RO*E1);

            if ( num <= 24 )
            {
                Gadgetron::maxAbsolute(mag, maxInten, ind);
            }
            else
            {
                hoNDArray<float> magPartial(RO, E1, 24, mag.get_data_ptr()+(num/2 - 12)*RO*E1);
                Gadgetron::maxAbsolute(magPartial, maxInten, ind);
            }
            if ( maxInten < FLT_EPSILON ) maxInten = 1.0f;

            if ( (maxInten<min_intensity_value_) || (maxInten>max_intensity_value_) )
            {
                GDEBUG_CONDITION_STREAM(verboseMode_, "Using the dynamic intensity scaling factor - may not have noise prewhitening performed ... ");
                scalingFactor_ = (float)(max_intensity_value_US_)/maxInten;
            }
            else
            {
                GDEBUG_CONDITION_STREAM(verboseMode_, "Using the fixed intensity scaling factor - must have noise prewhitening performed ... ");
                scalingFactor_ = SNR_NOISEFLOOR_SCALEFACTOR;

                while ( (maxInten*scalingFactor_ > max_intensity_value_) && (scalingFactor_>=2) )
                {
                    scalingFactor_ /= 2;
                }

                if (maxInten*scalingFactor_ > max_intensity_value_)
                {
                    GDEBUG_CONDITION_STREAM(verboseMode_, "The fixed intensity scaling factor leads to dynamic range overflow - switch to dyanmic intensity scaling ... ");
                    scalingFactor_ = (float)(max_intensity_value_)/maxInten;
                }

                use_constant_scalingFactor_ = true;
            }

            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_ : " << scalingFactor_);
            scal((float)scalingFactor_, mag);
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verboseMode_, "Using the fixed intensity scaling factor - scaling factor has been preset to be : " << scalingFactor_ << " ... ");
            scal((float)scalingFactor_, mag);
        }

        return true;
    }

    bool GtPlusReconGadget::
        generateKSpaceFilter(WorkOrderType& workOrder)
    {
        try
        {
            size_t RO = workOrder.data_.get_size(0);
            size_t E1 = workOrder.data_.get_size(1);
            size_t E2 = workOrder.data_.get_size(4);

            size_t RO_ref = workOrder.ref_.get_size(0);
            size_t E1_ref = workOrder.ref_.get_size(1);
            size_t E2_ref = workOrder.ref_.get_size(4);

            if ( workOrder.CalibMode_ == Gadgetron::ISMRMRD_interleaved )
            {
                RO_ref = RO;
                E1_ref = E1;
                E2_ref = E2;
            }

            // image data filter
            if ( RO>1 && filterRO_type_ != ISMRMRD_FILTER_NONE )
            {
                workOrder.filterRO_.create(RO);
                Gadgetron::generate_symmetric_filter(RO, workOrder.filterRO_, filterRO_type_, filterRO_sigma_, (size_t)std::ceil(filterRO_width_*RO));
                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterRO_, debugFolder_fullPath_+"filterRO"); }
            }

            if ( E1>1 && filterE1_type_ != ISMRMRD_FILTER_NONE )
            {
                // set up the filter lengh along E1
                if (max_E1_>0 && center_E1_>0)
                {
                    long long len = 2 * (center_E1_ - min_E1_);
                    if (len < 2 * (max_E1_ + 1 - center_E1_))
                    {
                        len = 2 * (max_E1_ + 1 - center_E1_);
                    }

                    if (E1>len)
                    {
                        workOrder.filterE1_.create(len);
                        Gadgetron::generate_symmetric_filter(len, workOrder.filterE1_, filterE1_type_, filterE1_sigma_, (size_t)std::ceil(filterE1_width_*len));

                        hoNDArray<ValueType> fE1;
                        Gadgetron::pad(E1, &workOrder.filterE1_, &fE1);
                        workOrder.filterE1_ = fE1;
                    }
                    else
                    {
                        workOrder.filterE1_.create(E1);
                        Gadgetron::generate_symmetric_filter(E1, workOrder.filterE1_, filterE1_type_, filterE1_sigma_, (size_t)std::ceil(filterE1_width_*E1));
                    }
                }
                else
                {
                    workOrder.filterE1_.create(E1);
                    Gadgetron::generate_symmetric_filter(E1, workOrder.filterE1_, filterE1_type_, filterE1_sigma_, (size_t)std::ceil(filterE1_width_*E1));
                }
            }

            if ( E2>1 && filterE2_type_ != ISMRMRD_FILTER_NONE )
            {
                if (max_E2_>0 && center_E2_>0)
                {
                    long long len = 2 * (center_E2_ - min_E2_);
                    if (len < 2 * (max_E2_ + 1 - center_E2_))
                    {
                        len = 2 * (max_E2_ + 1 - center_E2_);
                    }

                    if (E2>len)
                    {
                        workOrder.filterE2_.create(len);
                        Gadgetron::generate_symmetric_filter(len, workOrder.filterE2_, filterE2_type_, filterE2_sigma_, (size_t)std::ceil(filterE2_width_*len));

                        hoNDArray<ValueType> fE2;
                        Gadgetron::pad(E2, &workOrder.filterE2_, &fE2);
                        workOrder.filterE2_ = fE2;
                    }
                    else
                    {
                        workOrder.filterE2_.create(E2);
                        Gadgetron::generate_symmetric_filter(E2, workOrder.filterE2_, filterE2_type_, filterE2_sigma_, (size_t)std::ceil(filterE2_width_*E2));
                    }
                }
                else
                {
                    workOrder.filterE2_.create(E2);
                    Gadgetron::generate_symmetric_filter(E2, workOrder.filterE2_, filterE2_type_, filterE2_sigma_, (size_t)std::ceil(filterE2_width_*E2));
                }
            }

            // ref data filter
            if ( workOrder.ref_.get_number_of_elements() > 0 )
            {
                size_t startRO(0), endRO(0), startE1(0), endE1(0), startE2(0), endE2(0);
                if ( E2_ref == 1 )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlus_util_complex_.detectSampledRegion2D(workOrder.ref_, startRO, endRO, startE1, endE1));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlus_util_complex_.detectSampledRegion3D(workOrder.ref_, startRO, endRO, startE1, endE1, startE2, endE2));
                }

                if ( (workOrder.CalibMode_ == ISMRMRD_interleaved) || (workOrder.CalibMode_ == ISMRMRD_embedded) )
                {
                    // use the image data sample range
                    startRO = workOrder.start_RO_;
                    endRO = workOrder.end_RO_;
                }

                if ( RO_ref > 1 && filterRO_ref_type_ != ISMRMRD_FILTER_NONE )
                {
                    workOrder.filterRO_ref_.create(RO_ref);
                    Gadgetron::generate_symmetric_filter_ref(RO_ref, startRO, endRO, workOrder.filterRO_ref_);
                    if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterRO_ref_, debugFolder_fullPath_+"filterRO_ref"); }
                }

                if ( (workOrder.CalibMode_ == ISMRMRD_separate) || (workOrder.CalibMode_ == ISMRMRD_external) )
                {
                    if ( E1_ref > 1 && filterE1_ref_type_ != ISMRMRD_FILTER_NONE )
                    {
                        size_t len = endE1-startE1+1;
                        workOrder.filterE1_ref_.create(len);
                        Gadgetron::generate_symmetric_filter(len, workOrder.filterE1_ref_, filterE1_ref_type_, filterE1_ref_sigma_, (size_t)std::ceil(filterE1_ref_width_*len));
                        if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE1_ref_, debugFolder_fullPath_+"filterE1_ref"); }
                    }

                    if ( E2_ref > 1 && filterE2_ref_type_ != ISMRMRD_FILTER_NONE )
                    {
                        size_t len = endE2-startE2+1;
                        workOrder.filterE2_ref_.create(len);
                        Gadgetron::generate_symmetric_filter(len, workOrder.filterE2_ref_, filterE2_ref_type_, filterE2_ref_sigma_, (size_t)std::ceil(filterE2_ref_width_*len));
                        if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE2_ref_, debugFolder_fullPath_+"filterE2_ref"); }
                    }
                }
                else
                {
                    // this makes sure for interleaved and embedded, the kspace filter is applied at correct lines
                    if ( E1_ref > 1 && filterE1_ref_type_ != ISMRMRD_FILTER_NONE )
                    {
                        size_t len = E1_ref;
                        workOrder.filterE1_ref_.create(len);
                        Gadgetron::generate_symmetric_filter_ref(len, startE1, endE1, workOrder.filterE1_ref_);
                        if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE1_ref_, debugFolder_fullPath_+"filterE1_ref"); }
                    }

                    if ( E2_ref > 1 && filterE2_ref_type_ != ISMRMRD_FILTER_NONE )
                    {
                        size_t len = E2_ref;
                        workOrder.filterE2_ref_.create(len);
                        Gadgetron::generate_symmetric_filter_ref(len, startE2, endE2, workOrder.filterE2_ref_);
                        if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE2_ref_, debugFolder_fullPath_+"filterE2_ref"); }
                    }
                }
            }

            // partial fourier handling filter
            if ( RO>1 && workOrder.start_RO_>=0 && workOrder.end_RO_>0 )
            {
                workOrder.filterRO_partialfourier_.create(RO);
                Gadgetron::generate_asymmetric_filter(RO, workOrder.start_RO_, workOrder.end_RO_, workOrder.filterRO_partialfourier_, filterRO_pf_type_, (size_t)std::ceil(filterRO_pf_width_*RO), filterRO_pf_densityComp_);
                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterRO_partialfourier_, debugFolder_fullPath_+"filterRO_partialfourier"); }
            }

            if (E1>1 && workOrder.start_E1_ >= 0 && workOrder.end_E1_>0)
            {
                if (workOrder.start_E1_ == 0 || workOrder.end_E1_ == E1 - 1)
                {
                    workOrder.filterE1_partialfourier_.create(E1);
                    Gadgetron::generate_asymmetric_filter(E1, workOrder.start_E1_, workOrder.end_E1_, workOrder.filterE1_partialfourier_, filterE1_pf_type_, (size_t)std::ceil(filterE1_pf_width_*E1), filterE1_pf_densityComp_);
                }
                else
                {
                    size_t fil_len = workOrder.end_E1_ - workOrder.start_E1_ + 1;
                    size_t len_end = 2 * (workOrder.end_E1_ - E1 / 2 + 1);
                    size_t len_start = 2 * (E1 / 2 - workOrder.start_E1_ + 1);

                    if (len_end>len_start)
                    {
                        hoNDArray<ValueType> fil(len_end);
                        Gadgetron::generate_asymmetric_filter(len_end, len_end - fil_len, len_end - 1, fil, filterE1_pf_type_, (size_t)std::ceil(filterE1_pf_width_*len_end), filterE1_pf_densityComp_);
                        Gadgetron::pad(E1, &fil, &workOrder.filterE1_partialfourier_);
                    }
                    else
                    {
                        hoNDArray<ValueType> fil(len_start);
                        Gadgetron::generate_asymmetric_filter(len_start, 0, fil_len - 1, fil, filterE1_pf_type_, (size_t)std::ceil(filterE1_pf_width_*len_start), filterE1_pf_densityComp_);
                        Gadgetron::pad(E1, &fil, &workOrder.filterE1_partialfourier_);
                    }
                }

                if (!debugFolder_fullPath_.empty()) { gt_exporter_.exportArrayComplex(workOrder.filterE1_partialfourier_, debugFolder_fullPath_ + "filterE1_partialfourier"); }
            }

            if (E2>1 && workOrder.start_E2_ >= 0 && workOrder.end_E2_>0)
            {
                if (workOrder.start_E2_ == 0 || workOrder.end_E2_ == E2 - 1)
                {
                    workOrder.filterE2_partialfourier_.create(E2);
                    Gadgetron::generate_asymmetric_filter(E2, workOrder.start_E2_, workOrder.end_E2_, workOrder.filterE2_partialfourier_, filterE2_pf_type_, (size_t)std::ceil(filterE2_pf_width_*E2), filterE2_pf_densityComp_);
                }
                else
                {
                    size_t fil_len = workOrder.end_E2_ - workOrder.start_E2_ + 1;
                    size_t len_end = 2 * (workOrder.end_E2_ - E2 / 2 + 1);
                    size_t len_start = 2 * (E2 / 2 - workOrder.start_E2_ + 1);

                    if (len_end>len_start)
                    {
                        hoNDArray<ValueType> fil(len_end);
                        Gadgetron::generate_asymmetric_filter(len_end, len_end - fil_len, len_end - 1, fil, filterE2_pf_type_, (size_t)std::ceil(filterE2_pf_width_*len_end), filterE2_pf_densityComp_);
                        Gadgetron::pad(E2, &fil, &workOrder.filterE2_partialfourier_);
                    }
                    else
                    {
                        hoNDArray<ValueType> fil(len_start);
                        Gadgetron::generate_asymmetric_filter(len_start, 0, fil_len - 1, fil, filterE2_pf_type_, (size_t)std::ceil(filterE2_pf_width_*len_start), filterE2_pf_densityComp_);
                        Gadgetron::pad(E2, &fil, &workOrder.filterE2_partialfourier_);
                    }
                }

                if (!debugFolder_fullPath_.empty()) { gt_exporter_.exportArrayComplex(workOrder.filterE2_partialfourier_, debugFolder_fullPath_ + "filterE2_partialfourier"); }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors in GtPlusReconGadget::generateKSpaceFilter(...) ... ");
            return false;
        }

        return true;
    }

    bool GtPlusReconGadget::
        recomputeImageGeometry(GtPlusGadgetImageArray* images, GtPlusGadgetImageExt& imageHeader, size_t slc, size_t e2, size_t con, size_t phs, size_t rep, size_t set, size_t seg, size_t ave, size_t maxE2)
    {
        size_t E2 = images->matrix_size[4];

        // need to recompute image geometry
        // no need to consider RO and E1, because image position vector points to the image center

        if ( e2 >= E2 ) e2 = E2/2;

        size_t offsetCurr = images->get_offset(slc, e2, con, phs, rep, set, 0, ave);
        imageHeader = images->imageArray_[offsetCurr];

        // find the center partition
        if ( E2 > 1 )
        {
            long long midE2 = E2/2;
            size_t offset = images->get_offset(slc, midE2, con, phs, rep, set, 0, ave);

            while ( std::abs(imageHeader.slice_dir[0])<1e-6 && std::abs(imageHeader.slice_dir[1])<1e-6 && std::abs(imageHeader.slice_dir[2])<1e-6 )
            {
                imageHeader = images->imageArray_[offset];
                midE2++;
                offset = images->get_offset(slc, midE2, con, phs, rep, set, 0, ave);
            }

            // position vector for the center partition
            float posVec[3];
            posVec[0] = imageHeader.position[0];
            posVec[1] = imageHeader.position[1];
            posVec[2] = imageHeader.position[2];

            // slice direction
            float sliceVec[3];
            sliceVec[0] = imageHeader.slice_dir[0];
            sliceVec[1] = imageHeader.slice_dir[1];
            sliceVec[2] = imageHeader.slice_dir[2];

            midE2 = E2/2;

            // comput slice postion vector for this partition
            float posVecCurr[3];
            float e2_offset = (float)e2 - (float)midE2 + 0.5f;
            posVecCurr[0] = (float)(posVec[0] + aSpacing_[2] * sliceVec[0] * e2_offset);
            posVecCurr[1] = (float)(posVec[1] + aSpacing_[2] * sliceVec[1] * e2_offset);
            posVecCurr[2] = (float)(posVec[2] + aSpacing_[2] * sliceVec[2] * e2_offset);

            imageHeader.position[0] = posVecCurr[0];
            imageHeader.position[1] = posVecCurr[1];
            imageHeader.position[2] = posVecCurr[2];

            GDEBUG_CONDITION_STREAM(verboseMode_, "--> image position : [" << imageHeader.position[0] << " , " << imageHeader.position[1] << " , " << imageHeader.position[2] << "]");

            imageHeader.field_of_view[2] = (float)(aSpacing_[2]);

            imageHeader.user_int[0] = (int32_t)e2;
        }

        if ( imageHeader.measurement_uid == 0 )
        {
            GWARN_STREAM("imageHeader.measurement_uid == 0");
        }

        return true;
    }

    bool GtPlusReconGadget::
        sendOutRecon(GtPlusGadgetImageArray* images, const hoNDArray<ValueType>& res, int seriesNum, const std::vector<DimensionRecordType>& dimStartingIndexes, const std::string& prefix, const std::string& dataRole)
    {
        try
        {
            hoNDArray<real_value_type> timeStamp, physioTimeStamp;
            GADGET_CHECK_RETURN_FALSE( this->sendOutRecon(images, res, timeStamp, physioTimeStamp, seriesNum, dimStartingIndexes, prefix, dataRole) );
        }
        catch(...)
        {
            GERROR_STREAM("Errors in GtPlusReconGadget::sendOutRecon(complex float) ... ");
            return false;
        }

        return true;
    }

    bool GtPlusReconGadget::
        sendOutRecon(GtPlusGadgetImageArray* images, const hoNDArray<ValueType>& res, const hoNDArray<real_value_type>& timeStamp, const hoNDArray<real_value_type>& physioTimeStamp, 
        int seriesNum, const std::vector<DimensionRecordType>& dimStartingIndexes, const std::string& prefix, const std::string& dataRole)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dims = res.get_dimensions();
            size_t RO =  (*dims)[0];
            size_t E1 =  (*dims)[1];
            size_t CHA = (*dims)[2];
            size_t SLC = (*dims)[3];
            size_t E2 =  (*dims)[4];
            size_t CON = (*dims)[5];
            size_t PHS = (*dims)[6];
            size_t REP = (*dims)[7];
            size_t SET = (*dims)[8];
            size_t AVE = (*dims)[9];

            GDEBUG_CONDITION_STREAM(true, "sending out images, acquisition boundary [RO E1 CHA SLC E2 CON PHS REP SET AVE] = [" 
                << RO << " " << E1 << " " << CHA << " " 
                << SLC << " " << E2 << " " << CON << " " 
                << PHS << " " << REP << " " << SET << " " 
                << AVE << "] " );

            bool hasTimeStamp = false;
            if ( timeStamp.get_number_of_elements()>0 
                && timeStamp.get_size(9)==AVE 
                && timeStamp.get_size(8)==SET 
                && timeStamp.get_size(7)==REP 
                && timeStamp.get_size(6)==PHS 
                && timeStamp.get_size(5)==CON 
                && timeStamp.get_size(4)==E2 
                && timeStamp.get_size(3)==SLC )
            {
                hasTimeStamp = true;
            }

            bool hasPhysioTimeStamp = false;
            if ( physioTimeStamp.get_number_of_elements()>0 
                && physioTimeStamp.get_size(9)==AVE 
                && physioTimeStamp.get_size(8)==SET 
                && physioTimeStamp.get_size(7)==REP 
                && physioTimeStamp.get_size(6)==PHS 
                && physioTimeStamp.get_size(5)==CON 
                && physioTimeStamp.get_size(4)==E2 
                && physioTimeStamp.get_size(3)==SLC )
            {
                hasPhysioTimeStamp = true;
            }

            // info string for image, gfactor, snr map and std map
            std::ostringstream ostr_image;
            ostr_image << "x" << std::setprecision(4) << this->scalingFactor_;
            std::string imageInfo = ostr_image.str();

            std::ostringstream ostr_gfactor;
            ostr_gfactor << "x" << this->scalingFactor_gfactor_;
            std::string gfactorInfo = ostr_gfactor.str();

            std::ostringstream ostr_wrap_around_map;
            ostr_wrap_around_map << "x" << this->scalingFactor_wrap_around_map_;
            std::string wrapAroundMapInfo = ostr_wrap_around_map.str();

            std::ostringstream ostr_snr;
            ostr_snr << "x" << this->scalingFactor_snr_image_;
            std::string snrMapInfo = ostr_snr.str();

            std::ostringstream ostr_std;
            ostr_std << "x" << this->scalingFactor_std_map_;
            std::string stdMapInfo = ostr_std.str();

            // ------------------------------------------------------------- //

            std::vector<size_t> ind(10, 0);

            std::vector<size_t> dim2D(2);
            dim2D[0] = RO;
            dim2D[1] = E1;

            size_t set(0), rep(0), phs(0), con(0), e2(0), slc(0), cha(0), seg(0), ave(0);
            for ( ave=0; ave<AVE; ave++ )
            {
                for ( e2=0; e2<E2; e2++ )
                {
                    for ( slc=0; slc<SLC; slc++ )
                    {
                        for ( rep=0; rep<REP; rep++ )
                        {
                            for ( phs=0; phs<PHS; phs++ )
                            {
                                for ( set=0; set<SET; set++ )
                                {
                                    for ( con=0; con<CON; con++ )
                                    {
                                        GtPlusGadgetImageExt imageHeaderSent;

                                        GADGET_CHECK_RETURN_FALSE(recomputeImageGeometry(images, imageHeaderSent, slc, e2, con, phs, rep, set, 0, ave, E2));

                                        if ( imageHeaderSent.measurement_uid == 0 )
                                        {
                                            continue;
                                        }

                                        ind[0] = 0;
                                        ind[1] = 0;
                                        ind[2] = 0;
                                        ind[3] = slc;
                                        ind[4] = e2;
                                        ind[5] = con;
                                        ind[6] = phs;
                                        ind[7] = rep;
                                        ind[8] = set;
                                        ind[9] = ave;

                                        if ( hasTimeStamp )
                                        {
                                            if ( timeStamp(ind) > 0 )
                                            {
                                                imageHeaderSent.acquisition_time_stamp = (uint32_t)( (double)(timeStamp(ind)/timeStampResolution_) + 0.5 );
                                                GDEBUG_CONDITION_STREAM(verboseMode_, "Set acquisition time stamp : " << imageHeaderSent.acquisition_time_stamp);
                                            }
                                        }

                                        if ( hasPhysioTimeStamp )
                                        {
                                            if ( physioTimeStamp(ind) > 0 )
                                            {
                                                imageHeaderSent.physiology_time_stamp[0] = (uint32_t)( (double)(physioTimeStamp(ind)/timeStampResolution_) + 0.5 );
                                                GDEBUG_CONDITION_STREAM(verboseMode_, "Set physio time stamp : " << imageHeaderSent.physiology_time_stamp[0]);
                                            }
                                        }

                                        for ( cha=0; cha<CHA; cha++ )
                                        {
                                            ind[0] = 0;
                                            ind[1] = 0;
                                            ind[2] = cha;
                                            ind[3] = slc;
                                            ind[4] = e2;
                                            ind[5] = con;
                                            ind[6] = phs;
                                            ind[7] = rep;
                                            ind[8] = set;
                                            ind[9] = ave;

                                            hoNDArray<ValueType> currIm(dim2D, const_cast<ValueType*>(res.begin()+res.calculate_offset(ind)) );

                                            Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>();

                                            Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>* cm3 = new Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>();

                                            *(cm1->getObjectPtr()) = imageHeaderSent;

                                            cm1->getObjectPtr()->flags = 0;
                                            cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;

                                            // image number and image series
                                            cm1->getObjectPtr()->image_index = (uint16_t)computeSeriesImageNumber ( *(cm1->getObjectPtr()), CHA, cha, E2, e2);
                                            cm1->getObjectPtr()->image_series_index = seriesNum;
                                            // GDEBUG_CONDITION_STREAM(verboseMode_, "image number " << cm1->getObjectPtr()->image_index << "    image series " << cm1->getObjectPtr()->image_series_index << " ... ");

                                            // ----------------------------------------------------------
                                            // set the image attributes
                                            cm3->getObjectPtr()->set(GADGETRON_IMAGENUMBER, (long)cm1->getObjectPtr()->image_index);

                                            cm3->getObjectPtr()->set(GADGETRON_CHA,        (long)cha);
                                            cm3->getObjectPtr()->set(GADGETRON_SLC,        (long)cm1->getObjectPtr()->slice);
                                            cm3->getObjectPtr()->set(GADGETRON_E2,         (long)e2);
                                            cm3->getObjectPtr()->set(GADGETRON_CONTRAST,   (long)cm1->getObjectPtr()->contrast);
                                            cm3->getObjectPtr()->set(GADGETRON_PHASE,      (long)cm1->getObjectPtr()->phase);
                                            cm3->getObjectPtr()->set(GADGETRON_REP,        (long)cm1->getObjectPtr()->repetition);
                                            cm3->getObjectPtr()->set(GADGETRON_SET,        (long)cm1->getObjectPtr()->set);
                                            cm3->getObjectPtr()->set(GADGETRON_AVERAGE,    (long)cm1->getObjectPtr()->average);

                                            cm3->getObjectPtr()->set(GADGETRON_IMAGEPROCESSINGHISTORY, "GT");

                                            if ( dataRole == GADGETRON_IMAGE_REGULAR )
                                            {
                                                cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                                                cm3->getObjectPtr()->set(GADGETRON_IMAGECOMMENT, "GT");
                                                cm3->getObjectPtr()->append(GADGETRON_IMAGECOMMENT, imageInfo.c_str());

                                                cm3->getObjectPtr()->append(GADGETRON_SEQUENCEDESCRIPTION, "_GT");
                                                cm3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_REGULAR);
                                                cm3->getObjectPtr()->set(GADGETRON_IMAGE_SCALE_RATIO, (double)(this->scalingFactor_));
                                            }
                                            else if ( dataRole == GADGETRON_IMAGE_RETRO )
                                            {
                                                cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                                                cm3->getObjectPtr()->set(GADGETRON_IMAGECOMMENT, "GT");
                                                cm3->getObjectPtr()->append(GADGETRON_IMAGECOMMENT, "RETRO");
                                                cm3->getObjectPtr()->append(GADGETRON_IMAGECOMMENT, imageInfo.c_str());

                                                cm3->getObjectPtr()->set(GADGETRON_IMAGEPROCESSINGHISTORY, "RETRO");

                                                cm3->getObjectPtr()->set(GADGETRON_SEQUENCEDESCRIPTION, "_GT_RETRO");
                                                cm3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_RETRO);
                                                cm3->getObjectPtr()->set(GADGETRON_IMAGE_SCALE_RATIO, (double)(this->scalingFactor_));
                                            }
                                            else if ( dataRole == GADGETRON_IMAGE_PHASE )
                                            {
                                                cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_PHASE;

                                                cm3->getObjectPtr()->set(GADGETRON_IMAGECOMMENT, "PHS_GT");
                                                cm3->getObjectPtr()->set(GADGETRON_SEQUENCEDESCRIPTION, "PHS_GT");
                                                cm3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_PHASE);
                                                cm3->getObjectPtr()->set(GADGETRON_IMAGE_SCALE_RATIO, (double)(this->scalingFactor_));
                                            }
                                            else if ( dataRole == GADGETRON_IMAGE_GFACTOR )
                                            {
                                                cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                                                std::string comment = gfactorInfo;
                                                comment.append("_");
                                                comment.append("gfactor_GT");

                                                cm3->getObjectPtr()->set(GADGETRON_IMAGECOMMENT, comment.c_str());
                                                cm3->getObjectPtr()->set(GADGETRON_SEQUENCEDESCRIPTION, "_gfactor_GT");
                                                cm3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_GFACTOR);
                                                cm3->getObjectPtr()->set(GADGETRON_IMAGE_SCALE_RATIO, (double)(this->scalingFactor_gfactor_));
                                            }
                                            else if ( dataRole == GADGETRON_IMAGE_WRAPAROUNDMAP )
                                            {
                                                cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                                                std::string comment = wrapAroundMapInfo;
                                                comment.append("_");
                                                comment.append("WrapAround_Map_GT");

                                                cm3->getObjectPtr()->set(GADGETRON_IMAGECOMMENT, comment.c_str());
                                                cm3->getObjectPtr()->set(GADGETRON_SEQUENCEDESCRIPTION, "_WrapAround_Map_GT");
                                                cm3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_WRAPAROUNDMAP);
                                                cm3->getObjectPtr()->set(GADGETRON_IMAGE_SCALE_RATIO, (float)(this->scalingFactor_wrap_around_map_));
                                            }
                                            else if ( dataRole == GADGETRON_IMAGE_SNR_MAP )
                                            {
                                                cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                                                std::string comment = snrMapInfo;
                                                comment.append("_");
                                                comment.append("SNR_Map_GT");

                                                cm3->getObjectPtr()->set(GADGETRON_IMAGECOMMENT, comment.c_str());
                                                cm3->getObjectPtr()->set(GADGETRON_SEQUENCEDESCRIPTION, "_SNR_Map_GT");
                                                cm3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_SNR_MAP);
                                                cm3->getObjectPtr()->set(GADGETRON_IMAGE_SCALE_RATIO, (double)(this->scalingFactor_snr_image_));
                                            }
                                            else if ( dataRole == GADGETRON_IMAGE_STD_MAP )
                                            {
                                                cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                                                std::string comment = stdMapInfo;
                                                comment.append("_");
                                                comment.append("Std_Map_GT");

                                                cm3->getObjectPtr()->set(GADGETRON_IMAGECOMMENT, comment.c_str());
                                                cm3->getObjectPtr()->set(GADGETRON_SEQUENCEDESCRIPTION, "_Std_Map_GT");
                                                cm3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_STD_MAP);
                                                cm3->getObjectPtr()->set(GADGETRON_IMAGE_SCALE_RATIO, (double)(this->scalingFactor_std_map_));

                                                cm3->getObjectPtr()->set(GADGETRON_IMAGE_WINDOWCENTER, (long)(this->scalingFactor_std_map_));
                                                cm3->getObjectPtr()->set(GADGETRON_IMAGE_WINDOWWIDTH, (long)(2*this->scalingFactor_std_map_));
                                            }
                                            else if ( dataRole == GADGETRON_IMAGE_OTHER )
                                            {
                                                cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                                                cm3->getObjectPtr()->set(GADGETRON_IMAGECOMMENT, "GT");
                                                cm3->getObjectPtr()->set(GADGETRON_SEQUENCEDESCRIPTION, "_GT");
                                                cm3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_OTHER);
                                                cm3->getObjectPtr()->set(GADGETRON_IMAGE_SCALE_RATIO, (double)(this->scalingFactor_));
                                            }

                                            // ----------------------------------------------------------

                                            // set the time stamp
                                            // the time stamp of the first readout line in this 2D kspace is used

                                            Gadgetron::GadgetContainerMessage< Gadgetron::hoNDArray<ValueType> >* cm2 = new Gadgetron::GadgetContainerMessage< Gadgetron::hoNDArray<ValueType> >();
                                            cm1->cont(cm2);
                                            cm2->cont(cm3);

                                            std::vector<size_t> img_dims(2);
                                            img_dims[0] = RO;
                                            img_dims[1] = E1;

                                            //Fixing array dimensions (MSH)
                                            cm1->getObjectPtr()->matrix_size[0] = (uint16_t)RO;
                                            cm1->getObjectPtr()->matrix_size[1] = (uint16_t)E1;
                                            cm1->getObjectPtr()->matrix_size[2] = 1;
                                            cm1->getObjectPtr()->channels = 1;

                                            try
                                            {
                                                cm2->getObjectPtr()->create(&img_dims);
                                                Gadgetron::clear(cm2->getObjectPtr());
                                            }
                                            catch(...)
                                            {
                                                GDEBUG("Unable to allocate new image\n");
                                                cm1->release();
                                                return false;
                                            }

                                            memcpy(cm2->getObjectPtr()->begin(), currIm.begin(), sizeof(ValueType)*RO*E1);

                                            if ( !debugFolder2_fullPath_.empty() )
                                            {
                                                std::ostringstream ostr;
                                                ostr << prefix << "_" << cm1->getObjectPtr()->image_index;
                                                if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(*cm2->getObjectPtr(), debugFolder2_fullPath_+ostr.str()); }
                                            }

                                            GDEBUG_CONDITION_STREAM(verboseMode_, "sending out " << dataRole << " image [CHA SLC E2 CON PHS REP SET AVE] = [" 
                                                << cha << " " 
                                                << cm1->getObjectPtr()->slice << " " 
                                                << e2 << " " 
                                                << cm1->getObjectPtr()->contrast << " " 
                                                << cm1->getObjectPtr()->phase << " " 
                                                << cm1->getObjectPtr()->repetition << " " 
                                                << cm1->getObjectPtr()->set << " " 
                                                << cm1->getObjectPtr()->average << " " << "] " 
                                                << " -- Image number -- " << cm1->getObjectPtr()->image_index);

                                            // send out the images
                                            if (this->next()->putq(cm1) < 0) 
                                            {
                                                GERROR_STREAM("Put image to Q failed ... ");
                                                return false;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors in GtPlusReconGadget::sendOutRecon(complex float, time stamp) ... ");
            return false;
        }

        return true;
    }

    bool GtPlusReconGadget::sendOutRecon2D(GtPlusGadgetImageArray* images, const hoNDArray<ValueType>& res, int seriesNum, int imageNum)
    {
        try
        {
            // extract the magnitude
            hoNDArray<float> mag(res.get_dimensions());
            Gadgetron::abs(res, mag);
            GADGET_CHECK_RETURN_FALSE(scalingMagnitude(mag));
            GADGET_CHECK_RETURN_FALSE(sendOutRecon2D(images, mag, seriesNum, imageNum));
        }
        catch(...)
        {
            GERROR_STREAM("Exceptions happened in GtPlusReconGadget::sendOutRecon2D(...) ... ");
            return false;
        }

        return true;
    }

    bool GtPlusReconGadget::sendOutRecon2D(GtPlusGadgetImageArray* images, const hoNDArray<float>& res, int seriesNum, int imageNum)
    {
        try
        {
            Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new Gadgetron::GadgetContainerMessage<ISMRMRD::ImageHeader>();
            Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>* cm3 = new Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>();

            *(cm1->getObjectPtr()) = images->imageArray_[0];

            cm1->getObjectPtr()->flags = 0;
            cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;
            cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

            // image number and image series
            cm1->getObjectPtr()->image_index = imageNum;
            cm1->getObjectPtr()->image_series_index = seriesNum;

            Gadgetron::GadgetContainerMessage< Gadgetron::hoNDArray<float> >* cm2 = new Gadgetron::GadgetContainerMessage< Gadgetron::hoNDArray<float> >();
            cm1->cont(cm2);
            cm2->cont(cm3);

            std::vector<size_t> img_dims(2);
            img_dims[0] = res.get_size(0);
            img_dims[1] = res.get_size(1);

            // set the image attributes
            cm3->getObjectPtr()->set(GADGETRON_IMAGECOMMENT, "GT");
            cm3->getObjectPtr()->set(GADGETRON_SEQUENCEDESCRIPTION, "_GT");
            cm3->getObjectPtr()->set(GADGETRON_IMAGEPROCESSINGHISTORY, "GT");
            cm3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_REGULAR);

            cm3->getObjectPtr()->set(GADGETRON_CHA,        (long)0);
            cm3->getObjectPtr()->set(GADGETRON_SLC,        (long)cm1->getObjectPtr()->slice);
            cm3->getObjectPtr()->set(GADGETRON_E2,         (long)0);
            cm3->getObjectPtr()->set(GADGETRON_CONTRAST,   (long)cm1->getObjectPtr()->contrast);
            cm3->getObjectPtr()->set(GADGETRON_PHASE,      (long)cm1->getObjectPtr()->phase);
            cm3->getObjectPtr()->set(GADGETRON_REP,        (long)cm1->getObjectPtr()->repetition);
            cm3->getObjectPtr()->set(GADGETRON_SET,        (long)cm1->getObjectPtr()->set);
            cm3->getObjectPtr()->set(GADGETRON_AVERAGE,    (long)cm1->getObjectPtr()->average);

            cm3->getObjectPtr()->set(GADGETRON_IMAGE_SCALE_RATIO, (double)(this->scalingFactor_));

            //Fixing array dimensions (MSH)
            cm1->getObjectPtr()->matrix_size[0] = (uint16_t)res.get_size(0);
            cm1->getObjectPtr()->matrix_size[1] = (uint16_t)res.get_size(1);
            cm1->getObjectPtr()->matrix_size[2] = 1;
            cm1->getObjectPtr()->channels = 1;

            try
            {
                cm2->getObjectPtr()->create(&img_dims);
            }
            catch(...)
            {
                GDEBUG("Unable to allocate new image\n");
                cm1->release();
                return false;
            }

            memcpy(cm2->getObjectPtr()->begin(), res.begin(), sizeof(float)*res.get_size(0)*res.get_size(1));

            if ( !debugFolder2_fullPath_.empty() )
            {
                std::ostringstream ostr;
                ostr << "SentImage2D" << "_" << cm1->getObjectPtr()->image_index;
                if ( debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArray(*cm2->getObjectPtr(), debugFolder2_fullPath_+ostr.str()); }
            }

            // send out the images
            if (this->next()->putq(cm1) < 0) 
            {
                return false;
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors in GtPlusReconGadget::sendOutRecon2D(float) ... ");
            return false;
        }

        return true;
    }

    bool GtPlusReconGadget::computeSNRImage(const hoNDArray<ValueType>& res, const hoNDArray<ValueType>& gfactor, unsigned int startInd, bool withAcceleration, hoNDArray<ValueType>& snrImage, hoNDArray<ValueType>& stdMap)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dims = res.get_dimensions();
            size_t RO = (*dims)[0];
            size_t E1 = (*dims)[1];
            size_t CHA = (*dims)[2];
            size_t SLC = (*dims)[3];
            size_t E2 = (*dims)[4];
            size_t CON = (*dims)[5];
            size_t PHS = (*dims)[6];
            size_t REP = (*dims)[7];
            size_t SET = (*dims)[8];
            size_t AVE = (*dims)[9];

            snrImage = gfactor;

            if ( withAcceleration )
            {
                Gadgetron::addEpsilon(snrImage);
                Gadgetron::divide(res, snrImage, snrImage);
            }
            else
            {
                snrImage = res;
            }

            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(snrImage, debugFolder2_fullPath_+"snrImage"); }

            std::vector<size_t> dimStdMap(*dims);

            std::vector<size_t> ind(10, 0);
            size_t set(0), rep(0), phs(0), con(0), e2(0), slc(0), cha(0), seg(0), ave(0);

            if ( REP > startInd+2 )
            {
                dimStdMap[7] = 1;
                stdMap.create(dimStdMap);
                Gadgetron::clear(stdMap);

                size_t numOfIm = REP - startInd;

                hoNDArray<ValueType> repBuf(RO, E1, numOfIm);
                hoNDArray<real_value_type> repBufMag(RO, E1, numOfIm);
                hoNDArray<real_value_type> stdMap2D(RO, E1);

                for ( ave=0; ave<AVE; ave++ )
                {
                    for ( set=0; set<SET; set++ )
                    {
                        for ( phs=0; phs<PHS; phs++ )
                        {
                            for ( con=0; con<CON; con++ )
                            {
                                for ( e2=0; e2<E2; e2++ )
                                {
                                    for ( slc=0; slc<SLC; slc++ )
                                    {
                                        for ( cha=0; cha<CHA; cha++ )
                                        {
                                            Gadgetron::clear(repBuf);

                                            for ( rep=startInd; rep<REP; rep++ )
                                            {
                                                ind[2] = cha;
                                                ind[3] = slc;
                                                ind[4] = e2;
                                                ind[5] = con;
                                                ind[6] = phs;
                                                ind[7] = rep;
                                                ind[8] = set;
                                                ind[9] = ave;

                                                size_t offset = snrImage.calculate_offset(ind);

                                                memcpy(repBuf.begin()+(rep-startInd)*RO*E1, 
                                                    snrImage.begin()+offset, sizeof(ValueType)*RO*E1);
                                            }

                                            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(repBuf, debugFolder2_fullPath_+"repBuf"); }

                                            Gadgetron::abs(repBuf, repBufMag);
                                            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArray(repBufMag, debugFolder2_fullPath_+"repBufMag"); }

                                            // compute std
                                            GADGET_CHECK_RETURN_FALSE(Gadgetron::stdOver3rdDimension(repBufMag, stdMap2D, true));
                                            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArray(stdMap2D, debugFolder2_fullPath_+"stdMap2D"); }

                                            // copy it to the std map
                                            ind[2] = cha;
                                            ind[3] = slc;
                                            ind[4] = e2;
                                            ind[5] = con;
                                            ind[6] = phs;
                                            ind[7] = 0;
                                            ind[8] = set;
                                            ind[9] = ave;

                                            size_t offset = stdMap.calculate_offset(ind);
                                            hoNDArray<ValueType> stdMapCurr(RO, E1, stdMap.begin()+offset, false);
                                            stdMapCurr.copyFrom(stdMap2D);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else if ( PHS > startInd+2 )
            {
                dimStdMap[6] = 1;
                stdMap.create(dimStdMap);
                Gadgetron::clear(stdMap);

                size_t numOfIm = PHS - startInd;

                hoNDArray<ValueType> phsBuf(RO, E1, numOfIm);
                hoNDArray<real_value_type> phsBufMag(RO, E1, numOfIm);
                hoNDArray<real_value_type> stdMap2D(RO, E1);

                for ( ave=0; ave<AVE; ave++ )
                {
                    for ( set=0; set<SET; set++ )
                    {
                        for ( rep=0; rep<REP; rep++ )
                        {
                            for ( con=0; con<CON; con++ )
                            {
                                for ( e2=0; e2<E2; e2++ )
                                {
                                    for ( slc=0; slc<SLC; slc++ )
                                    {
                                        for ( cha=0; cha<CHA; cha++ )
                                        {
                                            Gadgetron::clear(phsBuf);

                                            for ( phs=startInd; phs<PHS; phs++ )
                                            {
                                                ind[2] = cha;
                                                ind[3] = slc;
                                                ind[4] = e2;
                                                ind[5] = con;
                                                ind[6] = phs;
                                                ind[7] = rep;
                                                ind[8] = set;
                                                ind[9] = ave;

                                                size_t offset = snrImage.calculate_offset(ind);

                                                memcpy(phsBuf.begin()+(phs-startInd)*RO*E1, 
                                                    snrImage.begin()+offset, sizeof(ValueType)*RO*E1);
                                            }

                                            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(phsBuf, debugFolder2_fullPath_+"phsBuf"); }

                                            Gadgetron::abs(phsBuf, phsBufMag);
                                            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArray(phsBufMag, debugFolder2_fullPath_+"phsBufMag"); }

                                            // compute std
                                            GADGET_CHECK_RETURN_FALSE(Gadgetron::stdOver3rdDimension(phsBufMag, stdMap2D, true));
                                            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArray(stdMap2D, debugFolder2_fullPath_+"stdMap2D"); }

                                            // copy it to the std map
                                            ind[2] = cha;
                                            ind[3] = slc;
                                            ind[4] = e2;
                                            ind[5] = con;
                                            ind[6] = 0;
                                            ind[7] = rep;
                                            ind[8] = set;
                                            ind[9] = ave;

                                            size_t offset = stdMap.calculate_offset(ind);
                                            hoNDArray<ValueType> stdMapCurr(RO, E1, stdMap.begin()+offset, false);
                                            stdMapCurr.copyFrom(stdMap2D);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors in GtPlusReconGadget::computeSNRImage(res, gfactor, snrImage, stdmap) ... ");
            return false;
        }

        return true;
    }

    int GtPlusReconGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GtPlusReconGadget - close(flags) : " << flags);

        if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

        if ( flags != 0 )
        {
            std::string procTime;
            gtPlus_util_.getCurrentMoment(procTime);

            GDEBUG_STREAM("* ============================================================================== *");
            GDEBUG_STREAM("---> MR recon phase, Currnt processing time : " << procTime << " <---");
            GDEBUG_STREAM("* ============================================================================== *");
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GtPlusReconGadget)

}
