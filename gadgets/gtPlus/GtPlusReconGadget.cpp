
#include "GtPlusReconGadget.h"
#include "GtPlusGadgetOpenMP.h"
#include "gadgetron_paths.h"
#include <iomanip>
#include "CloudBus.h"

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

    GtPlusReconGadget::GtPlusReconGadget() : mem_manager_(new Gadgetron::gtPlus::gtPlusMemoryManager(4, 640*1024*1024))
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

            min_intensity_value_ = this->get_int_value("min_intensity_value");
            GDEBUG_CONDITION_STREAM(verboseMode_, "min_intensity_value_ is " << min_intensity_value_);

            max_intensity_value_ = this->get_int_value("max_intensity_value");
            GDEBUG_CONDITION_STREAM(verboseMode_, "max_intensity_value_ is " << max_intensity_value_);

            scalingFactor_ = this->get_double_value("scalingFactor");
            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_ is " << scalingFactor_);

            scalingFactor_gfactor_ = this->get_double_value("scalingFactor_gfactor");
            if ( scalingFactor_gfactor_ == 0 ) scalingFactor_gfactor_ = 100;
            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_gfactor_ is " << scalingFactor_gfactor_);

            scalingFactor_wrap_around_map_ = this->get_double_value("scalingFactor_wrap_around_map");
            if ( scalingFactor_wrap_around_map_ == 0 ) scalingFactor_wrap_around_map_ = 1000;
            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_wrap_around_map_ is " << scalingFactor_wrap_around_map_);

            scalingFactor_snr_image_ = this->get_double_value("scalingFactor_snr_image");
            if ( scalingFactor_snr_image_ == 0 ) scalingFactor_snr_image_ = 10;
            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_snr_image_ is " << scalingFactor_snr_image_);

            scalingFactor_std_map_ = this->get_double_value("scalingFactor_std_map");
            if ( scalingFactor_std_map_ == 0 ) scalingFactor_std_map_ = 1000;
            GDEBUG_CONDITION_STREAM(verboseMode_, "scalingFactor_std_map_ is " << scalingFactor_std_map_);

            start_frame_for_std_map_ = this->get_int_value("start_frame_for_std_map");
            if ( start_frame_for_std_map_ == 0 ) start_frame_for_std_map_ = 5;
            GDEBUG_CONDITION_STREAM(verboseMode_, "start_frame_for_std_map_ is " << start_frame_for_std_map_);

            use_constant_scalingFactor_ = this->get_bool_value("use_constant_scalingFactor");
            GDEBUG_CONDITION_STREAM(verboseMode_, "use_constant_scalingFactor_ is " << use_constant_scalingFactor_);

            boost::shared_ptr<std::string> str = this->get_string_value("debugFolder");
            debugFolder_ = *str;
            GDEBUG_CONDITION_STREAM(verboseMode_, "debugFolder_ is " << debugFolder_);

            boost::shared_ptr<std::string> str2 = this->get_string_value("debugFolder2");
            debugFolder2_ = *str2;
            GDEBUG_CONDITION_STREAM(verboseMode_, "debugFolder2_ is " << debugFolder2_);

            timeStampResolution_ = (float)this->get_double_value("timeStampResolution");
            if ( timeStampResolution_ < FLT_EPSILON ) timeStampResolution_ = 0.0025f;
            GDEBUG_CONDITION_STREAM(verboseMode_, "timeStampResolution_ is " << timeStampResolution_);

            str = this->get_string_value("send_out_recon");
            if ( !str->empty() )
            {
                send_out_recon_ = this->get_bool_value("send_out_recon");
            }
            else
            {
                send_out_recon_ = true;
            }
            GDEBUG_CONDITION_STREAM(verboseMode_, "send_out_recon_ is " << send_out_recon_);

            str = this->get_string_value("send_out_recon_second");
            if ( !str->empty() )
            {
                send_out_recon_second_ = this->get_bool_value("send_out_recon_second");
            }
            else
            {
                send_out_recon_second_ = true;
            }
            GDEBUG_CONDITION_STREAM(verboseMode_, "send_out_recon_second_ is " << send_out_recon_second_);

            performTiming_ = this->get_bool_value("performTiming");
            GDEBUG_CONDITION_STREAM(verboseMode_, "performTiming_ is " << performTiming_);

            performTiming_ = this->get_bool_value("performTiming");
            GDEBUG_CONDITION_STREAM(verboseMode_, "performTiming_ is " << performTiming_);

            // kspace filter parameters
            str = this->get_string_value("filterRO");
            filterRO_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(*str);
            filterRO_sigma_ = this->get_double_value("filterRO_sigma");
            filterRO_width_ = this->get_double_value("filterRO_width");
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_type_ is " << *str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_sigma_ is " << filterRO_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_width_ is " << filterRO_width_);

            str = this->get_string_value("filterE1");
            filterE1_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(*str);
            filterE1_sigma_ = this->get_double_value("filterE1_sigma");
            filterE1_width_ = this->get_double_value("filterE1_width");
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_type_ is " << *str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_sigma_ is " << filterE1_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_width_ is " << filterE1_width_);

            str = this->get_string_value("filterE2");
            filterE2_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(*str);
            filterE2_sigma_ = this->get_double_value("filterE2_sigma");
            filterE2_width_ = this->get_double_value("filterE2_width");
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_type_ is " << *str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_sigma_ is " << filterE2_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_width_ is " << filterE2_width_);

            str = this->get_string_value("filterRefRO");
            filterRO_ref_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(*str);
            filterRO_ref_sigma_ = this->get_double_value("filterRefRO_sigma");
            filterRO_ref_width_ = this->get_double_value("filterRefRO_width");
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_ref_type_ is " << *str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_ref_sigma_ is " << filterRO_ref_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_ref_width_ is " << filterRO_ref_width_);

            str = this->get_string_value("filterRefE1");
            filterE1_ref_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(*str);
            filterE1_ref_sigma_ = this->get_double_value("filterRefE1_sigma");
            filterE1_ref_width_ = this->get_double_value("filterRefE1_width");
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_ref_type_ is " << *str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_ref_sigma_ is " << filterE1_ref_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_ref_width_ is " << filterE1_ref_width_);

            str = this->get_string_value("filterRefE2");
            filterE2_ref_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(*str);
            filterE2_ref_sigma_ = this->get_double_value("filterRefE2_sigma");
            filterE2_ref_width_ = this->get_double_value("filterRefE2_width");
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_ref_type_ is " << *str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_ref_sigma_ is " << filterE2_ref_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_ref_width_ is " << filterE2_ref_width_);

            str = this->get_string_value("filterPartialFourierRO");
            filterRO_pf_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(*str);
            filterRO_pf_sigma_ = this->get_double_value("filterPartialFourierRO_sigma");
            filterRO_pf_width_ = this->get_double_value("filterPartialFourierRO_width");
            filterRO_pf_densityComp_ = this->get_bool_value("filterPartialFourierRO_densityComp");
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_type_ is " << *str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_sigma_ is " << filterRO_pf_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_width_ is " << filterRO_pf_width_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_densityComp_ is " << filterRO_pf_densityComp_);

            str = this->get_string_value("filterPartialFourierE1");
            filterE1_pf_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(*str);
            filterE1_pf_sigma_ = this->get_double_value("filterPartialFourierE1_sigma");
            filterE1_pf_width_ = this->get_double_value("filterPartialFourierE1_width");
            filterE1_pf_densityComp_ = this->get_bool_value("filterPartialFourierE1_densityComp");
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_type_ is " << *str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_sigma_ is " << filterE1_pf_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_width_ is " << filterE1_pf_width_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_densityComp_ is " << filterE1_pf_densityComp_);

            str = this->get_string_value("filterPartialFourierE2");
            filterE2_pf_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(*str);
            filterE2_pf_sigma_ = this->get_double_value("filterPartialFourierE2_sigma");
            filterE2_pf_width_ = this->get_double_value("filterPartialFourierE2_width");
            filterE2_pf_densityComp_ = this->get_bool_value("filterPartialFourierE2_densityComp");
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_type_ is " << *str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_sigma_ is " << filterE2_pf_sigma_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_width_ is " << filterE2_pf_width_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_densityComp_ is " << filterE2_pf_densityComp_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            CloudComputing_ = this->get_bool_value("CloudComputing");
            CloudSize_ = (unsigned int)(this->get_int_value("CloudSize"));

            GDEBUG_CONDITION_STREAM(verboseMode_, "CloudComputing_ is " << CloudComputing_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "CloudSize_ is " << CloudSize_);

            str = this->get_string_value("cloudNodeFile");
            cloud_node_file_ = *str;
            GDEBUG_CONDITION_STREAM(verboseMode_, "cloud_node_file_ is " << cloud_node_file_);

            // read in the cloud information for every node
            gt_cloud_.resize(CloudSize_);

            for ( unsigned int ii=0; ii<CloudSize_; ii++ )
            {
                std::ostringstream ostreamstr1;
                ostreamstr1 << "CloudNode" << ii << "_IP" << std::ends;
                boost::shared_ptr<std::string> IP = this->get_string_value(ostreamstr1.str().c_str());
                gt_cloud_[ii].get<0>() = *IP;

                std::ostringstream ostreamstr2;
                ostreamstr2 << "CloudNode" << ii << "_Port" << std::ends;
                boost::shared_ptr<std::string> Port = this->get_string_value(ostreamstr2.str().c_str());
                gt_cloud_[ii].get<1>() = *Port;

                std::ostringstream ostreamstr3;
                ostreamstr3 << "CloudNode" << ii << "_XMLConfiguration" << std::ends;
                boost::shared_ptr<std::string> xmlName = this->get_string_value(ostreamstr3.str().c_str());
                gt_cloud_[ii].get<2>() = *xmlName;

                std::ostringstream ostreamstr4;
                ostreamstr4 << "CloudNode" << ii << "_ComputingPowerIndex" << std::ends;
                unsigned int computingPowerIndex = this->get_int_value(ostreamstr4.str().c_str());
                gt_cloud_[ii].get<3>() = computingPowerIndex;

                GDEBUG_CONDITION_STREAM(verboseMode_, "Cloud Node " << ii << " : " << gt_cloud_[ii]);
            }

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            thread_number_ratio_ = (float)this->get_double_value("thread_number_ratio");
            if ( thread_number_ratio_>1 || thread_number_ratio_<0 ) thread_number_ratio_ = 0;
            GDEBUG_CONDITION_STREAM(verboseMode_, "thread_number_ratio_ is " << thread_number_ratio_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "==================================================================");

            GDEBUG_CONDITION_STREAM(verboseMode_, "------> GtPlus recon parameters <------");

            workOrderPara_.upstream_coil_compression_ = this->get_bool_value("upstream_coil_compression");
            GDEBUG_CONDITION_STREAM(verboseMode_, "upstream_coil_compression_ is " << workOrderPara_.upstream_coil_compression_);

            workOrderPara_.upstream_coil_compression_thres_ = this->get_double_value("upstream_coil_compression_thres");
            GDEBUG_CONDITION_STREAM(verboseMode_, "upstream_coil_compression_thres_ is " << workOrderPara_.upstream_coil_compression_thres_);

            workOrderPara_.upstream_coil_compression_num_modesKept_ = this->get_int_value("upstream_coil_compression_num_modesKept");
            GDEBUG_CONDITION_STREAM(verboseMode_, "upstream_coil_compression_num_modesKept_ is " << workOrderPara_.upstream_coil_compression_num_modesKept_);

            workOrderPara_.downstream_coil_compression_ = this->get_bool_value("downstream_coil_compression");
            GDEBUG_CONDITION_STREAM(verboseMode_, "downstream_coil_compression_ is " << workOrderPara_.downstream_coil_compression_);

            workOrderPara_.coil_compression_thres_ = this->get_double_value("coil_compression_thres");

            if ( workOrderPara_.upstream_coil_compression_ && (workOrderPara_.coil_compression_thres_ > workOrderPara_.upstream_coil_compression_thres_) )
                workOrderPara_.coil_compression_thres_ = workOrderPara_.upstream_coil_compression_thres_;

            GDEBUG_CONDITION_STREAM(verboseMode_, "coil_compression_thres_ is " << workOrderPara_.coil_compression_thres_);

            workOrderPara_.coil_compression_num_modesKept_ = this->get_int_value("coil_compression_num_modesKept");

            if ( workOrderPara_.upstream_coil_compression_ && (workOrderPara_.coil_compression_num_modesKept_ > workOrderPara_.upstream_coil_compression_num_modesKept_) )
                workOrderPara_.coil_compression_num_modesKept_ = workOrderPara_.upstream_coil_compression_num_modesKept_;

            GDEBUG_CONDITION_STREAM(verboseMode_, "coil_compression_num_modesKept_ is " << workOrderPara_.coil_compression_num_modesKept_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            str = this->get_string_value("coil_map_algorithm");
            workOrderPara_.coil_map_algorithm_ = gtPlus_util_.getISMRMRDCoilMapAlgoFromName(*str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "coil_map_algorithm_ is " << *str);

            workOrderPara_.csm_kSize_ = (size_t)(this->get_int_value("csm_kSize"));
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_kSize_ is " << workOrderPara_.csm_kSize_);

            workOrderPara_.csm_powermethod_num_ = (size_t)(this->get_int_value("csm_powermethod_num"));
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_powermethod_num_ is " << workOrderPara_.csm_powermethod_num_);

            workOrderPara_.csm_true_3D_ = this->get_bool_value("csm_true_3D");
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_true_3D_ is " << workOrderPara_.csm_true_3D_);

            workOrderPara_.csm_iter_num_ = (size_t)(this->get_int_value("csm_iter_num"));
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_iter_num_ is " << workOrderPara_.csm_iter_num_);

            workOrderPara_.csm_iter_thres_ = this->get_double_value("csm_iter_thres");
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_iter_thres_ is " << workOrderPara_.csm_iter_thres_);

            workOrderPara_.csm_use_gpu_ = this->get_bool_value("csm_use_gpu");
            GDEBUG_CONDITION_STREAM(verboseMode_, "csm_use_gpu_ is " << workOrderPara_.csm_use_gpu_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            str = this->get_string_value("recon_algorithm");
            workOrderPara_.recon_algorithm_ = gtPlus_util_.getISMRMRDReconAlgoFromName(*str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "recon_algorithm_ is " << *str);

            workOrderPara_.recon_auto_parameters_ = this->get_bool_value("recon_auto_parameters");
            GDEBUG_CONDITION_STREAM(verboseMode_, "recon_auto_parameters_ is " << workOrderPara_.recon_auto_parameters_);

            workOrderPara_.gfactor_needed_ = this->get_bool_value("gfactor_needed");
            GDEBUG_CONDITION_STREAM(verboseMode_, "gfactor_needed_ is " << workOrderPara_.gfactor_needed_);

            workOrderPara_.wrap_around_map_needed_ = this->get_bool_value("wrap_around_map_needed");
            GDEBUG_CONDITION_STREAM(verboseMode_, "wrap_around_map_needed_ is " << workOrderPara_.wrap_around_map_needed_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.grappa_kSize_RO_ = (size_t)(this->get_int_value("grappa_kSize_RO"));
            workOrderPara_.grappa_kSize_E1_ = (size_t)(this->get_int_value("grappa_kSize_E1"));
            workOrderPara_.grappa_kSize_E2_ = (size_t)(this->get_int_value("grappa_kSize_E2"));
            workOrderPara_.grappa_reg_lamda_ = this->get_double_value("grappa_reg_lamda");
            workOrderPara_.grappa_calib_over_determine_ratio_ = this->get_double_value("grappa_calib_over_determine_ratio");
            workOrderPara_.grappa_use_gpu_ = this->get_bool_value("grappa_use_gpu");

            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_kSize_RO_ is " << workOrderPara_.grappa_kSize_RO_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_kSize_E1_ is " << workOrderPara_.grappa_kSize_E1_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_kSize_E2_ is " << workOrderPara_.grappa_kSize_E2_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_reg_lamda_ is " << workOrderPara_.grappa_reg_lamda_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_calib_over_determine_ratio_ is " << workOrderPara_.grappa_calib_over_determine_ratio_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "grappa_use_gpu_ is " << workOrderPara_.grappa_use_gpu_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.spirit_kSize_RO_ = (size_t)(this->get_int_value("spirit_kSize_RO"));
            if ( workOrderPara_.spirit_kSize_RO_ == 0 ) workOrderPara_.spirit_kSize_RO_ = 7;

            workOrderPara_.spirit_kSize_E1_ = (size_t)(this->get_int_value("spirit_kSize_E1"));
            if ( workOrderPara_.spirit_kSize_E1_ == 0 ) workOrderPara_.spirit_kSize_E1_ = 7;

            workOrderPara_.spirit_kSize_E2_ = (size_t)(this->get_int_value("spirit_kSize_E2"));
            if ( workOrderPara_.spirit_kSize_E2_ == 0 ) workOrderPara_.spirit_kSize_E2_ = 5;

            workOrderPara_.spirit_reg_lamda_ = this->get_double_value("spirit_reg_lamda");
            if ( workOrderPara_.spirit_reg_lamda_ < FLT_EPSILON ) workOrderPara_.spirit_reg_lamda_ = 0.005;

            workOrderPara_.spirit_use_gpu_ = this->get_bool_value("spirit_use_gpu");
            workOrderPara_.spirit_calib_over_determine_ratio_ = this->get_double_value("spirit_calib_over_determine_ratio");
            workOrderPara_.spirit_solve_symmetric_ = this->get_bool_value("spirit_solve_symmetric");

            workOrderPara_.spirit_iter_max_ = (size_t)(this->get_int_value("spirit_iter_max"));
            if ( workOrderPara_.spirit_iter_max_ == 0 ) workOrderPara_.spirit_iter_max_ = 100;

            workOrderPara_.spirit_iter_thres_ = this->get_double_value("spirit_iter_thres");
            if ( workOrderPara_.spirit_iter_thres_ < FLT_EPSILON ) workOrderPara_.spirit_iter_thres_ = 0.0015;

            workOrderPara_.spirit_print_iter_ = this->get_bool_value("spirit_print_iter");

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

            workOrderPara_.spirit_perform_linear_ = this->get_bool_value("spirit_perform_linear");
            workOrderPara_.spirit_perform_nonlinear_ = this->get_bool_value("spirit_perform_nonlinear");
            workOrderPara_.spirit_parallel_imaging_lamda_ = this->get_double_value("spirit_parallel_imaging_lamda");
            workOrderPara_.spirit_image_reg_lamda_ = this->get_double_value("spirit_image_reg_lamda");
            workOrderPara_.spirit_data_fidelity_lamda_ = this->get_double_value("spirit_data_fidelity_lamda");
            workOrderPara_.spirit_ncg_iter_max_ = (size_t)(this->get_int_value("spirit_ncg_iter_max"));
            workOrderPara_.spirit_ncg_iter_thres_ = this->get_double_value("spirit_ncg_iter_thres");
            workOrderPara_.spirit_ncg_print_iter_ = this->get_bool_value("spirit_ncg_print_iter");
            // spirit_ncg_scale_factor_ is computed from the data

            workOrderPara_.spirit_use_coil_sen_map_ = this->get_bool_value("spirit_use_coil_sen_map");
            workOrderPara_.spirit_use_moco_enhancement_ = this->get_bool_value("spirit_use_moco_enhancement");
            workOrderPara_.spirit_recon_moco_images_ = this->get_bool_value("spirit_recon_moco_images");
            workOrderPara_.spirit_RO_enhancement_ratio_ = this->get_double_value("spirit_RO_enhancement_ratio");
            workOrderPara_.spirit_E1_enhancement_ratio_ = this->get_double_value("spirit_E1_enhancement_ratio");
            workOrderPara_.spirit_E2_enhancement_ratio_ = this->get_double_value("spirit_E2_enhancement_ratio");
            workOrderPara_.spirit_temporal_enhancement_ratio_ = this->get_double_value("spirit_temporal_enhancement_ratio");
            workOrderPara_.spirit_2D_scale_per_chunk_ = this->get_bool_value("spirit_2D_scale_per_chunk");
            workOrderPara_.spirit_3D_scale_per_chunk_ = this->get_bool_value("spirit_3D_scale_per_chunk");

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

            str = this->get_string_value("retro_gated_interp_method");
            workOrderPara_.retro_gated_interp_method_ = gtPlus_util_.getISMRMRDRetroGatingInterpFromName(*str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "retro_gated_interp_method_ is " << *str);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            workOrderPara_.job_split_by_S_ = this->get_bool_value("job_split_by_S");
            workOrderPara_.job_num_of_N_ = (size_t)(this->get_int_value("job_num_of_N"));
            workOrderPara_.job_max_Megabytes_ = (size_t)(this->get_int_value("job_max_Megabytes"));
            workOrderPara_.job_overlap_ = (size_t)(this->get_int_value("job_overlap"));
            workOrderPara_.job_perform_on_control_node_ = this->get_bool_value("job_perform_on_control_node");

            GDEBUG_CONDITION_STREAM(verboseMode_, "job_split_by_S_ is " << workOrderPara_.job_split_by_S_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "job_num_of_N_ is " << workOrderPara_.job_num_of_N_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "job_max_Megabytes_ is " << workOrderPara_.job_max_Megabytes_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "job_overlap_ is " << workOrderPara_.job_overlap_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "job_perform_on_control_node_ is " << workOrderPara_.job_perform_on_control_node_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            str = this->get_string_value("partialFourier_algo");
            workOrderPara_.partialFourier_algo_ = gtPlus_util_.getISMRMRDPartialFourierReconAlgoFromName(*str);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_algo_ is " << *str);

            workOrderPara_.partialFourier_homodyne_iters_ = (size_t)(this->get_int_value("partialFourier_homodyne_iters"));
            workOrderPara_.partialFourier_homodyne_thres_ = this->get_double_value("partialFourier_homodyne_thres");
            workOrderPara_.partialFourier_homodyne_densityComp_ = this->get_bool_value("partialFourier_homodyne_densityComp");

            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_homodyne_iters_ is " << workOrderPara_.partialFourier_homodyne_iters_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_homodyne_thres_ is " << workOrderPara_.partialFourier_homodyne_thres_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_homodyne_densityComp_ is " << workOrderPara_.partialFourier_homodyne_densityComp_);

            workOrderPara_.partialFourier_POCS_iters_ = (size_t)(this->get_int_value("partialFourier_POCS_iters"));
            workOrderPara_.partialFourier_POCS_thres_ = this->get_double_value("partialFourier_POCS_thres");
            workOrderPara_.partialFourier_POCS_transitBand_ = (size_t)(this->get_int_value("partialFourier_POCS_transitBand"));
            workOrderPara_.partialFourier_POCS_transitBand_E2_ = (size_t)(this->get_int_value("partialFourier_POCS_transitBand_E2"));

            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_POCS_iters_ is " << workOrderPara_.partialFourier_POCS_iters_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_POCS_thres_ is " << workOrderPara_.partialFourier_POCS_thres_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_POCS_transitBand_ is " << workOrderPara_.partialFourier_POCS_transitBand_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_POCS_transitBand_ is " << workOrderPara_.partialFourier_POCS_transitBand_E2_);

            workOrderPara_.partialFourier_FengHuang_kSize_RO_ = (size_t)(this->get_int_value("partialFourier_FengHuang_kSize_RO"));
            workOrderPara_.partialFourier_FengHuang_kSize_E1_ = (size_t)(this->get_int_value("partialFourier_FengHuang_kSize_E1"));
            workOrderPara_.partialFourier_FengHuang_kSize_E2_ = (size_t)(this->get_int_value("partialFourier_FengHuang_kSize_E2"));
            workOrderPara_.partialFourier_FengHuang_thresReg_ = this->get_double_value("partialFourier_FengHuang_thresReg");
            workOrderPara_.partialFourier_FengHuang_sameKernel_allN_ = this->get_bool_value("partialFourier_FengHuang_sameKernel_allN");
            workOrderPara_.partialFourier_FengHuang_transitBand_ = (size_t)(this->get_int_value("partialFourier_FengHuang_transitBand"));
            workOrderPara_.partialFourier_FengHuang_transitBand_E2_ = (size_t)(this->get_int_value("partialFourier_FengHuang_transitBand_E2"));

            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_kSize_RO_ is " << workOrderPara_.partialFourier_FengHuang_kSize_RO_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_kSize_E1_ is " << workOrderPara_.partialFourier_FengHuang_kSize_E1_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_kSize_E2_ is " << workOrderPara_.partialFourier_FengHuang_kSize_E2_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_thresReg_ is " << workOrderPara_.partialFourier_FengHuang_thresReg_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_sameKernel_allN_ is " << workOrderPara_.partialFourier_FengHuang_sameKernel_allN_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_transitBand_ is " << workOrderPara_.partialFourier_FengHuang_transitBand_);
            GDEBUG_CONDITION_STREAM(verboseMode_, "partialFourier_FengHuang_transitBand_E2_ is " << workOrderPara_.partialFourier_FengHuang_transitBand_E2_);

            GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

            recon_kspace_needed_ = this->get_bool_value("recon_kspace_needed");
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

        bool using_cloudbus = this->get_bool_value("using_cloudbus");
        bool has_cloud_node_xml_configuration = this->get_string_value("CloudNodeXMLConfiguration")->size();

        if (using_cloudbus && has_cloud_node_xml_configuration) {
            std::vector<GadgetronNodeInfo> nodes;
            CloudBus::instance()->get_node_info(nodes);
            gtCloud.resize(nodes.size());

            unsigned int n;
            for ( n=0; n<nodes.size(); n++ )
            {
                std::stringstream ss;
                gtCloud[n].get<0>() = nodes[n].address;
                ss << nodes[n].port;
                gtCloud[n].get<1>() = ss.str();
                gtCloud[n].get<2>() = *this->get_string_value("CloudNodeXMLConfiguration");
                gtCloud[n].get<3>() = nodes[n].compute_capability;

                GDEBUG_CONDITION_STREAM(verboseMode_, "Gadget Node " << n << " : " << gt_cloud_[n]);
            }

            return true; //We will leave the function here

        }

        std::string nodeFileName = get_gadgetron_home();
        nodeFileName.append("/config/gtCloud/");
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

        verboseMode_ = this->get_bool_value("verboseMode");

        // read parameters from xml
        image_series_ = this->get_int_value("image_series");

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
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilter(RO, workOrder.start_RO_, workOrder.end_RO_, workOrder.filterRO_, filterRO_type_, filterRO_sigma_, (size_t)std::ceil(filterRO_width_*RO)));
                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterRO_, debugFolder_fullPath_+"filterRO"); }
            }

            if ( E1>1 && filterE1_type_ != ISMRMRD_FILTER_NONE )
            {
                workOrder.filterE1_.create(E1);
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilter(E1, workOrder.start_E1_, workOrder.end_E1_, workOrder.filterE1_, filterE1_type_, filterE1_sigma_, (size_t)std::ceil(filterE1_width_*E1)));
                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE1_, debugFolder_fullPath_+"filterE1"); }
            }

            if ( E2>1 && filterE2_type_ != ISMRMRD_FILTER_NONE )
            {
                workOrder.filterE2_.create(E2);
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilter(E2, workOrder.start_E2_, workOrder.end_E2_, workOrder.filterE2_, filterE2_type_, filterE2_sigma_, (size_t)std::ceil(filterE2_width_*E2)));
                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE2_, debugFolder_fullPath_+"filterE2"); }
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
                    startRO = workOrder.start_RO_; if ( startRO < 0 ) startRO=0;
                    endRO = workOrder.end_RO_; if ( endRO < 0 ) endRO = RO_ref-1;
                }

                if ( RO_ref > 1 && filterRO_ref_type_ != ISMRMRD_FILTER_NONE )
                {
                    workOrder.filterRO_ref_.create(RO_ref);
                    GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilterForRef(RO_ref, startRO, endRO, workOrder.filterRO_ref_, filterRO_ref_type_, filterRO_ref_sigma_, (size_t)std::ceil(filterRO_ref_width_*RO_ref)));
                    if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterRO_ref_, debugFolder_fullPath_+"filterRO_ref"); }
                }

                if ( (workOrder.CalibMode_ == ISMRMRD_separate) || (workOrder.CalibMode_ == ISMRMRD_external) )
                {
                    if ( E1_ref > 1 && filterE1_ref_type_ != ISMRMRD_FILTER_NONE )
                    {
                        size_t len = endE1-startE1+1;
                        workOrder.filterE1_ref_.create(len);
                        GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilter(len, 0, len-1, workOrder.filterE1_ref_, filterE1_ref_type_, filterE1_ref_sigma_, (size_t)std::ceil(filterE1_ref_width_*len)));
                        if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE1_ref_, debugFolder_fullPath_+"filterE1_ref"); }
                    }

                    if ( E2_ref > 1 && filterE2_ref_type_ != ISMRMRD_FILTER_NONE )
                    {
                        size_t len = endE2-startE2+1;
                        workOrder.filterE2_ref_.create(len);
                        GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilter(len, 0, len-1, workOrder.filterE2_ref_, filterE2_ref_type_, filterE2_ref_sigma_, (size_t)std::ceil(filterE2_ref_width_*len)));
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
                        GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilterForRef(len, startE1, endE1, workOrder.filterE1_ref_, filterE1_ref_type_, filterE1_ref_sigma_, (size_t)std::ceil(filterE1_ref_width_*len)));
                        if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE1_ref_, debugFolder_fullPath_+"filterE1_ref"); }
                    }

                    if ( E2_ref > 1 && filterE2_ref_type_ != ISMRMRD_FILTER_NONE )
                    {
                        size_t len = E2_ref;
                        workOrder.filterE2_ref_.create(len);
                        GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilterForRef(len, startE2, endE2, workOrder.filterE2_ref_, filterE2_ref_type_, filterE2_ref_sigma_, (size_t)std::ceil(filterE2_ref_width_*len)));
                        if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE2_ref_, debugFolder_fullPath_+"filterE2_ref"); }
                    }
                }
            }

            // partial fourier handling filter
            if ( RO>1 && workOrder.start_RO_>=0 && workOrder.end_RO_>0 )
            {
                workOrder.filterRO_partialfourier_.create(RO);
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateAsymmetricFilter(RO, workOrder.start_RO_, workOrder.end_RO_, workOrder.filterRO_partialfourier_, filterRO_pf_type_, (size_t)std::ceil(filterRO_pf_width_*RO), filterRO_pf_densityComp_));
                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterRO_partialfourier_, debugFolder_fullPath_+"filterRO_partialfourier"); }
            }

            if ( E1>1 && workOrder.start_E1_>=0 && workOrder.end_E1_>0 )
            {
                workOrder.filterE1_partialfourier_.create(E1);
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateAsymmetricFilter(E1, workOrder.start_E1_, workOrder.end_E1_, workOrder.filterE1_partialfourier_, filterE1_pf_type_, (size_t)std::ceil(filterE1_pf_width_*E1), filterE1_pf_densityComp_));
                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE1_partialfourier_, debugFolder_fullPath_+"filterE1_partialfourier"); }
            }

            if ( E2>1 && workOrder.start_E2_>=0 && workOrder.end_E2_>0 )
            {
                workOrder.filterE2_partialfourier_.create(E2);
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateAsymmetricFilter(E2, workOrder.start_E2_, workOrder.end_E2_, workOrder.filterE2_partialfourier_, filterE2_pf_type_, (size_t)std::ceil(filterE2_pf_width_*E2), filterE2_pf_densityComp_));
                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE2_partialfourier_, debugFolder_fullPath_+"filterE2_partialfourier"); }
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
            size_t midE2 = E2/2;
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
            posVecCurr[0] = (float)(posVec[0] + aSpacing_[2]*sliceVec[0]*(e2-midE2+0.5f));
            posVecCurr[1] = (float)(posVec[1] + aSpacing_[2]*sliceVec[1]*(e2-midE2+0.5f));
            posVecCurr[2] = (float)(posVec[2] + aSpacing_[2]*sliceVec[2]*(e2-midE2+0.5f));

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
