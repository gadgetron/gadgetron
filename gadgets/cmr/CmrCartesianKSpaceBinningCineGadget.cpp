
#include "CmrCartesianKSpaceBinningCineGadget.h"

namespace Gadgetron {

    CmrCartesianKSpaceBinningCineGadget::CmrCartesianKSpaceBinningCineGadget() : BaseClass()
    {
        send_out_multiple_series_by_slice_ = false;
    }

    CmrCartesianKSpaceBinningCineGadget::~CmrCartesianKSpaceBinningCineGadget()
    {
    }

    int CmrCartesianKSpaceBinningCineGadget::process_config(const mrd::Header& header)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(header) == GADGET_OK, GADGET_FAIL);

        this->send_out_multiple_series_by_slice_ = this->send_out_multiple_series_by_slice.value();

        // -------------------------------------------------

        auto& h = header;

        if (h.user_parameters)
        {
            for (std::vector<mrd::UserParameterLongType>::const_iterator i = h.user_parameters->user_parameter_long.begin(); i != h.user_parameters->user_parameter_long.end(); ++i)
            {
                if (std::strcmp(i->name.c_str(), "MultiSeriesForSlices") == 0)
                {
                    GDEBUG_STREAM("Found from protocol, MultiSeriesForSlices is defined ... ");
                    this->send_out_multiple_series_by_slice_ = true;
                    GDEBUG_STREAM("Reset, send_out_multiple_series_by_slice_ is " << send_out_multiple_series_by_slice_);
                }
            }
        }

        // -------------------------------------------------

        binning_reconer_.debug_folder_                                   = this->debug_folder_full_path_;
        binning_reconer_.perform_timing_                                 = this->perform_timing.value();
        binning_reconer_.verbose_                                        = this->verbose.value();

        binning_reconer_.use_multiple_channel_recon_                     = this->use_multiple_channel_recon.value();
        binning_reconer_.use_paralell_imaging_binning_recon_             = true;
        binning_reconer_.use_nonlinear_binning_recon_                    = this->use_nonlinear_binning_recon.value();

        binning_reconer_.estimate_respiratory_navigator_                 = true;
        binning_reconer_.respiratory_navigator_moco_reg_strength_        = this->respiratory_navigator_moco_reg_strength.value();
        binning_reconer_.respiratory_navigator_moco_iters_               = this->respiratory_navigator_moco_iters.value();

        binning_reconer_.time_tick_                                      = this->time_tick.value();
        binning_reconer_.trigger_time_index_                             = 0;
        binning_reconer_.arrhythmia_rejector_factor_                     = this->arrhythmia_rejector_factor.value();

        binning_reconer_.grappa_kSize_RO_                                = this->grappa_kSize_RO.value();
        binning_reconer_.grappa_kSize_E1_                                = this->grappa_kSize_E1.value();
        binning_reconer_.grappa_reg_lamda_                               = this->grappa_reg_lamda.value();
        binning_reconer_.downstream_coil_compression_num_modesKept_      = this->downstream_coil_compression_num_modesKept.value();
        binning_reconer_.downstream_coil_compression_thres_              = this->downstream_coil_compression_thres.value();

        binning_reconer_.kspace_binning_interpolate_heart_beat_images_   = this->kspace_binning_interpolate_heart_beat_images.value();
        binning_reconer_.kspace_binning_navigator_acceptance_window_     = this->kspace_binning_navigator_acceptance_window.value();

        binning_reconer_.kspace_binning_moco_reg_strength_               = this->kspace_binning_moco_reg_strength.value();
        binning_reconer_.kspace_binning_moco_iters_                      = this->kspace_binning_moco_iters.value();

        binning_reconer_.kspace_binning_max_temporal_window_             = this->kspace_binning_max_temporal_window.value();
        binning_reconer_.kspace_binning_minimal_cardiac_phase_width_     = this->kspace_binning_minimal_cardiac_phase_width.value();
        binning_reconer_.kspace_binning_kSize_RO_                        = this->kspace_binning_kSize_RO.value();
        binning_reconer_.kspace_binning_kSize_E1_                        = this->kspace_binning_kSize_E1.value();
        binning_reconer_.kspace_binning_reg_lamda_                       = this->kspace_binning_reg_lamda.value();
        binning_reconer_.kspace_binning_linear_iter_max_                 = this->kspace_binning_linear_iter_max.value();
        binning_reconer_.kspace_binning_linear_iter_thres_               = this->kspace_binning_linear_iter_thres.value();
        binning_reconer_.kspace_binning_nonlinear_iter_max_              = this->kspace_binning_nonlinear_iter_max.value();
        binning_reconer_.kspace_binning_nonlinear_iter_thres_            = this->kspace_binning_nonlinear_iter_thres.value();
        binning_reconer_.kspace_binning_nonlinear_data_fidelity_lamda_   = this->kspace_binning_nonlinear_data_fidelity_lamda.value();
        binning_reconer_.kspace_binning_nonlinear_image_reg_lamda_       = this->kspace_binning_nonlinear_image_reg_lamda.value();
        binning_reconer_.kspace_binning_nonlinear_reg_N_weighting_ratio_ = this->kspace_binning_nonlinear_reg_N_weighting_ratio.value();
        binning_reconer_.kspace_binning_nonlinear_reg_use_coil_sen_map_  = this->kspace_binning_nonlinear_reg_use_coil_sen_map.value();
        binning_reconer_.kspace_binning_nonlinear_reg_with_approx_coeff_ = this->kspace_binning_nonlinear_reg_with_approx_coeff.value();
        binning_reconer_.kspace_binning_nonlinear_reg_wav_name_          = this->kspace_binning_nonlinear_reg_wav_name.value();

        return GADGET_OK;
    }

    int CmrCartesianKSpaceBinningCineGadget::process(Gadgetron::GadgetContainerMessage< mrd::ReconData >* m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("CmrCartesianKSpaceBinningCineGadget::process"); }

        process_called_times_++;

        mrd::ReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->buffers.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->buffers.size() << " instead of " << num_encoding_spaces_);
        }

        std::vector<unsigned int> processed_slices = kspace_binning_processed_slices.value();
        if(processed_slices.size()>0)
        {
            size_t ii;
            for (ii=0; ii<recon_bit_->buffers[0].data.headers.get_number_of_elements(); ii++)
            {
                if( recon_bit_->buffers[0].data.headers(ii).acquisition_time_stamp>0 ) break;
            }

            size_t curr_slc = recon_bit_->buffers[0].data.headers(ii).idx.slice.value_or(0);

            GDEBUG_STREAM("Incoming slice : " << curr_slc);

            bool do_processing = false;
            for (size_t k=0; k<processed_slices.size(); k++)
            {
                if(curr_slc==processed_slices[k])
                {
                    do_processing = true;
                    break;
                }
            }

            if(!do_processing)
            {
                GDEBUG_STREAM("Ignore incoming slice : " << curr_slc);
                m1->release();
                return GADGET_OK;
            }
        }

        // for every encoding space
        for (size_t e = 0; e < recon_bit_->buffers.size(); e++)
        {
            std::stringstream os;
            os << "_encoding_" << e;

            GDEBUG_CONDITION_STREAM(verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(verbose.value(), "======================================================================");

            // ---------------------------------------------------------------
            // export incoming data

            /*if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(recon_bit_->buffers[e].data.data, debug_folder_full_path_ + "data" + os.str());
            }

            if (!debug_folder_full_path_.empty() && recon_bit_->buffers[e].data.trajectory_)
            {
                if (recon_bit_->buffers[e].ref_->trajectory_->get_number_of_elements() > 0)
                {
                    gt_exporter_.export_array(*(recon_bit_->buffers[e].data.trajectory_), debug_folder_full_path_ + "data_traj" + os.str());
                }
            }*/

            // ---------------------------------------------------------------

            if (recon_bit_->buffers[e].data.data.get_number_of_elements() > 0)
            {
                /*if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array_complex(recon_bit_->buffers[e].data.data, debug_folder_full_path_ + "data_before_unwrapping" + os.str());
                }

                if (!debug_folder_full_path_.empty() && recon_bit_->buffers[e].data.trajectory_)
                {
                    if (recon_bit_->buffers[e].data.trajectory_->get_number_of_elements() > 0)
                    {
                        gt_exporter_.export_array(*(recon_bit_->buffers[e].data.trajectory_), debug_folder_full_path_ + "data_before_unwrapping_traj" + os.str());
                    }
                }*/

                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("CmrCartesianKSpaceBinningCineGadget::perform_binning"); }
                this->perform_binning(recon_bit_->buffers[e], e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("CmrCartesianKSpaceBinningCineGadget::compute_image_header, raw images"); }
                this->compute_image_header(recon_bit_->buffers[e], res_raw_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                this->set_time_stamps(res_raw_, acq_time_raw_, cpt_time_raw_);

                // ---------------------------------------------------------------

                if (!debug_folder_full_path_.empty())
                {
                    this->gt_exporter_.export_array_complex(res_raw_.data, debug_folder_full_path_ + "recon_res_raw" + os.str());
                }

                if(this->send_out_raw.value())
                {
                    if (perform_timing.value()) { gt_timer_.start("CmrCartesianKSpaceBinningCineGadget::send_out_image_array, raw"); }
                    this->send_out_image_array(res_raw_, e, image_series.value() + ((int)e + 1), GADGETRON_IMAGE_REGULAR);
                    if (perform_timing.value()) { gt_timer_.stop(); }
                }

                // ---------------------------------------------------------------
                this->create_binning_image_headers_from_raw();
                this->set_time_stamps(res_binning_, acq_time_binning_, cpt_time_binning_);

                if (!debug_folder_full_path_.empty())
                {
                    this->gt_exporter_.export_array_complex(res_binning_.data, debug_folder_full_path_ + "recon_res_binning" + os.str());
                }

                if (perform_timing.value()) { gt_timer_.start("CmrCartesianKSpaceBinningCineGadget::send_out_image_array, binning"); }
                this->send_out_image_array(res_binning_, e, image_series.value() + (int)e + 2, GADGETRON_IMAGE_RETRO);
                if (perform_timing.value()) { gt_timer_.stop(); }
            }
        }

        m1->release();

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    void CmrCartesianKSpaceBinningCineGadget::perform_binning(mrd::ReconAssembly& recon_bit, size_t encoding)
    {
        try
        {
            size_t RO  = recon_bit.data.data.get_size(0);
            size_t E1  = recon_bit.data.data.get_size(1);
            size_t E2  = recon_bit.data.data.get_size(2);
            size_t CHA = recon_bit.data.data.get_size(3);
            size_t N   = recon_bit.data.data.get_size(4);
            size_t S   = recon_bit.data.data.get_size(5);
            size_t SLC = recon_bit.data.data.get_size(6);

            size_t binned_N = this->number_of_output_phases.value();

            GADGET_CHECK_THROW(E2==1);
            GADGET_CHECK_THROW(N>binned_N);

            Gadgetron::GadgetronTimer timer(false);

            res_raw_.data.create(RO, E1, E2, 1, N, S, SLC);
            acq_time_raw_.create(N, S, SLC);
            cpt_time_raw_.create(N, S, SLC);

            res_binning_.data.create(RO, E1, E2, 1, binned_N, S, SLC);
            acq_time_binning_.create(binned_N, S, SLC);
            cpt_time_binning_.create(binned_N, S, SLC);

            size_t n, s, slc;
            for (slc=0; slc<SLC; slc++)
            {
                std::stringstream os;

                size_t ind = 0;
                while (recon_bit.data.headers[ind].measurement_uid==0 && ind< recon_bit.data.headers.get_number_of_elements())
                {
                    ind++;
                }

                size_t curr_slc = recon_bit.data.headers[ind].idx.slice.value_or(0);

                os << "_encoding_" << encoding << "_SLC_" << slc << "_SLCOrder_" << curr_slc;

                std::string suffix = os.str();

                GDEBUG_STREAM("Processing binning on SLC : " << slc << " - " << curr_slc << " , encoding space : " << encoding << " " << suffix);

                // set up the binning object
                binning_reconer_.binning_obj_.data_.create(RO, E1, CHA, N, S, recon_bit.data.data.begin()+slc*RO*E1*CHA*N*S);
                binning_reconer_.binning_obj_.sampling_ = recon_bit.data.sampling;
                binning_reconer_.binning_obj_.headers_.create(E1, N, S, recon_bit.data.headers.begin()+slc*E1*N*S);

                binning_reconer_.binning_obj_.output_N_ = binned_N;
                binning_reconer_.binning_obj_.accel_factor_E1_ = acceFactorE1_[encoding];
                binning_reconer_.binning_obj_.random_sampling_ = (calib_mode_[encoding]!=mrd::CalibrationMode::kEmbedded
                                                                && calib_mode_[encoding]!=mrd::CalibrationMode::kInterleaved
                                                                && calib_mode_[encoding]!=mrd::CalibrationMode::kSeparate
                                                                && calib_mode_[encoding]!=mrd::CalibrationMode::kNoacceleration);

                binning_reconer_.suffix_ = suffix;

                // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(binning_reconer_.binning_obj_.data, debug_folder_full_path_ + "binning_obj_data" + os.str()); }

                // compute the binning
                if (perform_timing.value()) { timer.start("compute binning ... "); }
                try
                {
                    binning_reconer_.process_binning_recon();
                }
                catch(...)
                {
                    GERROR_STREAM("Exceptions happened in binning_reconer_.process_binning_recon() for slice " << slc);
                    continue;
                }
                if (perform_timing.value()) { timer.stop(); }

                if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(binning_reconer_.binning_obj_.complex_image_raw_, debug_folder_full_path_ + "binning_obj_complex_image_raw" + os.str()); }
                if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(binning_reconer_.binning_obj_.complex_image_binning_, debug_folder_full_path_ + "binning_obj_complex_image_binning" + os.str()); }

                // get the binnig results
                memcpy(this->res_raw_.data.begin()+slc*RO*E1*N*S*SLC,
                        binning_reconer_.binning_obj_.complex_image_raw_.begin(),
                        binning_reconer_.binning_obj_.complex_image_raw_.get_number_of_bytes());

                memcpy(this->res_binning_.data.begin()+slc*RO*E1*binned_N*S*SLC,
                        binning_reconer_.binning_obj_.complex_image_binning_.begin(),
                        binning_reconer_.binning_obj_.complex_image_binning_.get_number_of_bytes());

                for (s=0; s<S; s++)
                {
                    for (n=0; n<N; n++)
                    {
                        acq_time_raw_(n, s, slc) = binning_reconer_.binning_obj_.phs_time_stamp_(n, s);
                        cpt_time_raw_(n, s, slc) = binning_reconer_.binning_obj_.phs_cpt_time_stamp_(n, s);
                    }

                    for (n=0; n<binned_N; n++)
                    {
                        acq_time_binning_(n, s, slc) = binning_reconer_.binning_obj_.phs_time_stamp_(n, s);
                        cpt_time_binning_(n, s, slc) = binning_reconer_.binning_obj_.mean_RR_ * binning_reconer_.binning_obj_.desired_cpt_[n];
                    }
                }
            }

            std::stringstream os;
            os << "_encoding_" << encoding;

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(res_raw_.data, debug_folder_full_path_ + "binning_complex_image_raw" + os.str()); }
            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(res_binning_.data, debug_folder_full_path_ + "binning_complex_image_binning_" + os.str()); }

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array(acq_time_raw_, debug_folder_full_path_ + "binning_acq_time_raw" + os.str()); }
            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array(cpt_time_raw_, debug_folder_full_path_ + "binning_cpt_time_raw" + os.str()); }

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array(acq_time_binning_, debug_folder_full_path_ + "binning_acq_time_binning" + os.str()); }
            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array(cpt_time_binning_, debug_folder_full_path_ + "binning_cpt_time_binning" + os.str()); }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in CmrCartesianKSpaceBinningCineGadget::perform_binning(...) ... ");
        }
    }

    void CmrCartesianKSpaceBinningCineGadget::create_binning_image_headers_from_raw()
    {
        try
        {
            size_t N = res_raw_.headers.get_size(0);
            size_t S = res_raw_.headers.get_size(1);
            size_t SLC = res_raw_.headers.get_size(2);

            size_t binned_N = this->number_of_output_phases.value();

            res_binning_.headers.create(binned_N, S, SLC);
            res_binning_.meta.create(binned_N, S, SLC);

            size_t n, s, slc;
            for (slc=0; slc<SLC; slc++)
            {
                for (s=0; s<S; s++)
                {
                    for (n=0; n<binned_N; n++)
                    {
                        res_binning_.headers(n, s, slc) = res_raw_.headers(n, s, slc);
                        res_binning_.meta[n + s*binned_N + slc*binned_N*S] = res_raw_.meta[n + s*N + slc*N*S];
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in CmrCartesianKSpaceBinningCineGadget::create_binning_image_headers_from_raw(...) ... ");
        }
    }

    void CmrCartesianKSpaceBinningCineGadget::set_time_stamps(mrd::ImageArray& res, hoNDArray< float >& acq_time, hoNDArray< float >& cpt_time)
    {
        try
        {
            size_t N = res.headers.get_size(0);
            size_t S = res.headers.get_size(1);
            size_t SLC = res.headers.get_size(2);

            size_t n, s, slc;
            for (slc=0; slc<SLC; slc++)
            {
                for (s=0; s<S; s++)
                {
                    for (n=0; n<N; n++)
                    {
                        res.headers(n, s, slc).acquisition_time_stamp = (uint32_t)(acq_time(n, s, slc)/ this->time_tick.value() + 0.5);
                        res.headers(n, s, slc).physiology_time_stamp[0] = (uint32_t)(cpt_time(n, s, slc)/ this->time_tick.value() + 0.5);
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in CmrCartesianKSpaceBinningCineGadget::set_time_stamps(...) ... ");
        }
    }

    int CmrCartesianKSpaceBinningCineGadget::prep_image_header_send_out(mrd::ImageArray& res, size_t n, size_t s, size_t slc, size_t encoding, int series_num, const std::string& data_role)
    {
        try
        {
            ImageArraySendMixin::prep_image_header_send_out(res, n, s, slc, encoding, series_num, data_role);

            size_t RO = res.data.get_size(0);
            size_t E1 = res.data.get_size(1);
            size_t E2 = res.data.get_size(2);
            size_t CHA = res.data.get_size(3);
            size_t N = res.data.get_size(4);
            size_t S = res.data.get_size(5);
            size_t SLC = res.data.get_size(6);

            if(this->send_out_multiple_series_by_slice_)
            {
                res.headers(n, s, slc).image_series_index = res.headers(n, s, slc).image_series_index.value_or(0) + 100 * res.headers(n, s, slc).slice.value_or(0);

                size_t offset = n + s*N + slc*N*S;

                std::ostringstream ostr;
                ostr << "_SLC" << res.headers(n, s, slc).slice.value_or(0)+1;
                res.meta[offset][GADGETRON_SEQUENCEDESCRIPTION].push_back(ostr.str().c_str());
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors in GenericReconGadget::prep_image_header_send_out(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(CmrCartesianKSpaceBinningCineGadget)
}
