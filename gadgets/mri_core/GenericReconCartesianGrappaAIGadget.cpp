
#include "GenericReconCartesianGrappaAIGadget.h"
#include "mri_core_grappa.h"
#include "hoNDArray_reductions.h"

namespace Gadgetron {

    GenericReconCartesianGrappaAIGadget::GenericReconCartesianGrappaAIGadget() : BaseClass()
    {
        
        
    }

    GenericReconCartesianGrappaAIGadget::~GenericReconCartesianGrappaAIGadget()
    {
    }

    int GenericReconCartesianGrappaAIGadget::process_config(ACE_Message_Block *mb)
    {
        boost::filesystem::path gadgetron_python_path
            = this->context.paths.gadgetron_home / "share" / "gadgetron" / "python";
        GDEBUG_STREAM("PYTHON_PATH " << gadgetron_python_path.string());
        Gadgetron::initialize_python();
        Gadgetron::add_python_path(gadgetron_python_path.generic_string());
        this->gt_home_ = gadgetron_python_path.generic_string();

        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        // -------------------------------------------------

        models_.resize(num_encoding_spaces_);
        kernels_.resize(num_encoding_spaces_);
        recon_res_grappa_ai_.resize(num_encoding_spaces_);

        return GADGET_OK;
    }

    int GenericReconCartesianGrappaAIGadget::process(Gadgetron::GadgetContainerMessage<IsmrmrdReconData> *m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("GenericReconCartesianGrappaAIGadget::process"); }

        process_called_times_++;

        IsmrmrdReconData *recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size()
                << " instead of "
                << num_encoding_spaces_);
        }

        GadgetContainerMessage< std::vector<ISMRMRD::Waveform> > * wav = AsContainerMessage< std::vector<ISMRMRD::Waveform>  >(m1->cont());
        if (wav)
        {
            if (verbose.value())
            {
                GDEBUG_STREAM("Incoming recon_bit with " << wav->getObjectPtr()->size() << " wave form samples ");
            }
        }

        // for every encoding space
        for (size_t e = 0; e < recon_bit_->rbit_.size(); e++)
        {
            std::stringstream os;
            os << "_encoding_" << e << "_" << process_called_times_;

            GDEBUG_CONDITION_STREAM(verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(verbose.value(), "======================================================================");

            // ---------------------------------------------------------------
            // export incoming data

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_, debug_folder_full_path_ + "data" + os.str()); }
            if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].data_.trajectory_)
            {
                if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0) { gt_exporter_.export_array(*(recon_bit_->rbit_[e].data_.trajectory_), debug_folder_full_path_ + "data_traj" + os.str()); }
            }

            // ---------------------------------------------------------------

            if (recon_bit_->rbit_[e].ref_)
            {
                if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_bit_->rbit_[e].ref_->data_, debug_folder_full_path_ + "ref" + os.str()); }

                if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].ref_->trajectory_)
                {
                    if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0) { gt_exporter_.export_array(*(recon_bit_->rbit_[e].ref_->trajectory_), debug_folder_full_path_ + "ref_traj" + os.str()); }
                }

                // ---------------------------------------------------------------

                // after this step, the recon_obj_[e].ref_calib_ and recon_obj_[e].ref_coil_map_ are set

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianGrappaAIGadget::make_ref_coil_map"); }
                this->make_ref_coil_map(*recon_bit_->rbit_[e].ref_, *recon_bit_->rbit_[e].data_.data_.get_dimensions(), recon_obj_[e].ref_calib_, recon_obj_[e].ref_coil_map_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ----------------------------------------------------------
                // export prepared ref for calibration and coil map
                if (!debug_folder_full_path_.empty()) { this->gt_exporter_.export_array_complex(recon_obj_[e].ref_calib_, debug_folder_full_path_ + "ref_calib" + os.str()); }
                if (!debug_folder_full_path_.empty()) { this->gt_exporter_.export_array_complex(recon_obj_[e].ref_coil_map_, debug_folder_full_path_ + "ref_coil_map" + os.str()); }

                // ---------------------------------------------------------------
                // after this step, the recon_obj_[e].ref_calib_dst_ and recon_obj_[e].ref_coil_map_ are modified
                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianGrappaAIGadget::prepare_down_stream_coil_compression_ref_data"); }
                this->prepare_down_stream_coil_compression_ref_data(recon_obj_[e].ref_calib_, recon_obj_[e].ref_coil_map_, recon_obj_[e].ref_calib_dst_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                if (!debug_folder_full_path_.empty()) { this->gt_exporter_.export_array_complex(recon_obj_[e].ref_calib_dst_, debug_folder_full_path_ + "ref_calib_dst" + os.str()); }
                if (!debug_folder_full_path_.empty()) { this->gt_exporter_.export_array_complex(recon_obj_[e].ref_coil_map_, debug_folder_full_path_ + "ref_coil_map_dst" + os.str()); }

                // ---------------------------------------------------------------

                // after this step, coil map is computed and stored in recon_obj_[e].coil_map_
                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianGrappaAIGadget::perform_coil_map_estimation"); }
                this->perform_coil_map_estimation(recon_obj_[e].ref_coil_map_, recon_obj_[e].coil_map_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                // after this step, recon_obj_[e].kernel_, recon_obj_[e].kernelIm_, recon_obj_[e].unmixing_coeff_ are filled
                // gfactor is computed too
                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianGrappaAIGadget::perform_calib"); }
                this->perform_calib(recon_bit_->rbit_[e], recon_obj_[e], e);
                if (perform_timing.value()) { gt_timer_.stop(); }
                // ---------------------------------------------------------------

                recon_bit_->rbit_[e].ref_->clear();
                recon_bit_->rbit_[e].ref_ = Core::none;
            }

            if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0)
            {
                if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_, debug_folder_full_path_ + "data_before_unwrapping" + os.str()); }
                if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].data_.trajectory_)
                {
                    if (recon_bit_->rbit_[e].data_.trajectory_->get_number_of_elements() > 0) { gt_exporter_.export_array(*(recon_bit_->rbit_[e].data_.trajectory_), debug_folder_full_path_ + "data_before_unwrapping_traj" + os.str()); }
                }

                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianGrappaAIGadget::perform_unwrapping"); }
                this->perform_unwrapping(recon_bit_->rbit_[e], recon_obj_[e], e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianGrappaAIGadget::compute_image_header"); }
                this->compute_image_header(recon_bit_->rbit_[e], recon_obj_[e].recon_res_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // set up recon_res_grappa_ai_
                this->recon_res_grappa_ai_[e].headers_ = recon_obj_[e].recon_res_.headers_;
                this->recon_res_grappa_ai_[e].meta_ = recon_obj_[e].recon_res_.meta_;

                for (size_t m = 0; m < this->recon_res_grappa_ai_[e].meta_.size(); m++)
                {
                    this->recon_res_grappa_ai_[e].meta_[m].append(GADGETRON_IMAGECOMMENT, GADGETRON_AI);
                    this->recon_res_grappa_ai_[e].meta_[m].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_AI);
                }

                // ---------------------------------------------------------------
                // pass down waveform
                if (wav)
                {
                    recon_obj_[e].recon_res_.waveform_ = *wav->getObjectPtr();
                    this->recon_res_grappa_ai_[e].waveform_ = *wav->getObjectPtr();
                }

                recon_obj_[e].recon_res_.acq_headers_ = recon_bit_->rbit_[e].data_.headers_;
                this->recon_res_grappa_ai_[e].acq_headers_ = recon_bit_->rbit_[e].data_.headers_;

                // ---------------------------------------------------------------

                if (!debug_folder_full_path_.empty()) { this->gt_exporter_.export_array_complex(recon_obj_[e].recon_res_.data_, debug_folder_full_path_ + "recon_res" + os.str()); }
                if (!debug_folder_full_path_.empty()) { this->gt_exporter_.export_array_complex(this->recon_res_grappa_ai_[e].data_, debug_folder_full_path_ + "recon_res_grappa_ai" + os.str()); }

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianGrappaAIGadget::send_out_image_array"); }
                this->send_out_image_array(this->recon_res_grappa_ai_[e], e, image_series.value() + ((int)e + 2), GADGETRON_IMAGE_REGULAR);
                if (perform_timing.value()) { gt_timer_.stop(); }
            }

            recon_obj_[e].recon_res_.data_.clear();
            recon_obj_[e].gfactor_.clear();
            recon_obj_[e].recon_res_.headers_.clear();
            recon_obj_[e].recon_res_.meta_.clear();

            recon_res_grappa_ai_[e].data_.clear();
            recon_res_grappa_ai_[e].headers_.clear();
            recon_res_grappa_ai_[e].meta_.clear();
        }

        m1->release();

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    void GenericReconCartesianGrappaAIGadget::perform_calib(IsmrmrdReconBit &recon_bit, ReconObjType &recon_obj, size_t e)
    {
        size_t RO = recon_bit.data_.data_.get_size(0);
        size_t E1 = recon_bit.data_.data_.get_size(1);
        size_t E2 = recon_bit.data_.data_.get_size(2);

        hoNDArray<T> &src = recon_obj.ref_calib_;
        hoNDArray<T> &dst = recon_obj.ref_calib_dst_;

        size_t ref_RO = src.get_size(0);
        size_t ref_E1 = src.get_size(1);
        size_t ref_E2 = src.get_size(2);
        size_t srcCHA = src.get_size(3);
        size_t ref_N = src.get_size(4);
        size_t ref_S = src.get_size(5);
        size_t ref_SLC = src.get_size(6);

        models_[e].resize(ref_N*ref_S*ref_SLC);
        kernels_[e].resize(ref_N*ref_S*ref_SLC);

        size_t dstCHA = dst.get_size(3);

        recon_obj.unmixing_coeff_.create(RO, E1, E2, srcCHA, ref_N, ref_S, ref_SLC);
        recon_obj.gfactor_.create(RO, E1, E2, 1, ref_N, ref_S, ref_SLC);

        Gadgetron::clear(recon_obj.unmixing_coeff_);
        Gadgetron::clear(recon_obj.gfactor_);

        if (acceFactorE1_[e] <= 1 && acceFactorE2_[e] <= 1)
        {
            Gadgetron::conjugate(recon_obj.coil_map_, recon_obj.unmixing_coeff_);
        }
        else
        {
            // allocate buffer for kernels
            size_t kRO = grappa_kSize_RO.value();
            size_t kNE1 = grappa_kSize_E1.value();
            size_t kNE2 = grappa_kSize_E2.value();

            bool fitItself = this->downstream_coil_compression.value();

            size_t convKRO(1), convKE1(1), convKE2(1);

            if (E2 > 1)
            {
                std::vector<int> kE1, oE1;
                std::vector<int> kE2, oE2;
                grappa3d_kerPattern(kE1, oE1, kE2, oE2, convKRO, convKE1, convKE2, (size_t)acceFactorE1_[e],
                    (size_t)acceFactorE2_[e], kRO, kNE1, kNE2, fitItself);
            }
            else
            {
                std::vector<int> kE1, oE1;
                Gadgetron::grappa2d_kerPattern(kE1, oE1, convKRO, convKE1, (size_t)acceFactorE1_[e], kRO, kNE1,
                    fitItself);
                recon_obj.kernelIm_.create(RO, E1, 1, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);
            }

            recon_obj.kernel_.create(convKRO, convKE1, convKE2, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);

            Gadgetron::clear(recon_obj.kernel_);
            Gadgetron::clear(recon_obj.kernelIm_);

            // declare the grappa training
            PythonFunction<boost::python::object> create_grappa_ai_model("grappa_ai", "create_grappa_ai_model");

            // loop over every [N S SLC] for calibration
            long long num = ref_N * ref_S * ref_SLC;

            long long ii;
            for (ii = 0; ii < num; ii++)
            {
                size_t slc = ii / (ref_N * ref_S);
                size_t s = (ii - slc * ref_N * ref_S) / (ref_N);
                size_t n = ii - slc * ref_N * ref_S - s * ref_N;

                std::stringstream os;
                os << "n" << n << "_s" << s << "_slc" << slc << "_encoding_" << e;
                std::string suffix = os.str();

                T*pSrc = &(src(0, 0, 0, 0, n, s, slc));
                hoNDArray<T> ref_src(ref_RO, ref_E1, ref_E2, srcCHA, pSrc);

                T*pDst = &(dst(0, 0, 0, 0, n, s, slc));
                hoNDArray<T> ref_dst(ref_RO, ref_E1, ref_E2, dstCHA, pDst);

                // ---------------------------------------------------------------------

                if (E2 > 1)
                {
                    hoNDArray<T> ker(convKRO, convKE1, convKE2, srcCHA, dstCHA,
                        &(recon_obj.kernel_(0, 0, 0, 0, 0, n, s, slc)));

                    if (fitItself)
                    {
                        Gadgetron::grappa3d_calib_convolution_kernel(ref_src, ref_dst, (size_t)acceFactorE1_[e],
                            (size_t)acceFactorE2_[e], grappa_reg_lamda.value(),
                            grappa_calib_over_determine_ratio.value(), kRO, kNE1,
                            kNE2, ker);
                    }
                    else
                    {
                        Gadgetron::grappa3d_calib_convolution_kernel(ref_src, ref_src, (size_t)acceFactorE1_[e],
                            (size_t)acceFactorE2_[e], grappa_reg_lamda.value(),
                            grappa_calib_over_determine_ratio.value(), kRO, kNE1,
                            kNE2, ker);
                    }

                    if (!debug_folder_full_path_.empty())
                    {
                        gt_exporter_.export_array_complex(ker, debug_folder_full_path_ + "convKer3D_" + suffix);
                    }

                    hoNDArray<T> coilMap(RO, E1, E2, dstCHA,
                        &(recon_obj.coil_map_(0, 0, 0, 0, n, s, slc)));
                    hoNDArray<T> unmixC(RO, E1, E2, srcCHA,
                        &(recon_obj.unmixing_coeff_(0, 0, 0, 0, n, s, slc)));
                    hoNDArray<float> gFactor(RO, E1, E2, 1, &(recon_obj.gfactor_(0, 0, 0, 0, n, s, slc)));
                    Gadgetron::grappa3d_unmixing_coeff(ker, coilMap, (size_t)acceFactorE1_[e],
                        (size_t)acceFactorE2_[e], unmixC, gFactor);

                    if (!debug_folder_full_path_.empty())
                    {
                        gt_exporter_.export_array_complex(unmixC, debug_folder_full_path_ + "unmixC_3D_" + suffix);
                    }

                    if (!debug_folder_full_path_.empty())
                    {
                        gt_exporter_.export_array(gFactor, debug_folder_full_path_ + "gFactor_3D_" + suffix);
                    }
                }
                else
                {
                    hoNDArray<T> acsSrc(ref_RO, ref_E1, srcCHA,
                        const_cast<T*>(ref_src.begin()));
                    hoNDArray<T> acsDst(ref_RO, ref_E1, dstCHA,
                        const_cast<T*>(ref_dst.begin()));

                    hoNDArray<T> convKer(convKRO, convKE1, srcCHA, dstCHA,
                        &(recon_obj.kernel_(0, 0, 0, 0, 0, n, s, slc)));
                    hoNDArray<T> kIm(RO, E1, srcCHA, dstCHA,
                        &(recon_obj.kernelIm_(0, 0, 0, 0, 0, n, s, slc)));

                    size_t startRO = 0;
                    size_t endRO = acsSrc.get_size(0) - 1;
                    size_t startE1 = 0;
                    size_t endE1 = acsSrc.get_size(1) - 1;

                    std::vector<int> kE1, oE1;

                    // ----------------------------------------
                    // typical grappa calibration
                    // compute ker, A and B
                    // ----------------------------------------

                    size_t convkRO, convkE1;

                    Gadgetron::grappa2d_kerPattern(kE1, oE1, convkRO, convkE1, (size_t)acceFactorE1_[e], kRO, kNE1, fitItself);

                    hoNDArray< T> ker;
                    hoNDArray< T> A, B;

                    if (fitItself)
                    {
                        Gadgetron::grappa2d_prepare_calib(acsSrc, acsDst, kRO, kE1, oE1, startRO, endRO, startE1, endE1, A, B);
                    }
                    else
                    {
                        Gadgetron::grappa2d_prepare_calib(acsSrc, acsSrc, kRO, kE1, oE1, startRO, endRO, startE1, endE1, A, B);
                    }

                    double thres = grappa_reg_lamda.value();
                    Gadgetron::grappa2d_perform_calib(A, B, kRO, kE1, oE1, thres, ker);

                    kernels_[e][ii] = ker;

                    if (!debug_folder_full_path_.empty())
                    {
                        gt_exporter_.export_array_complex(ker, debug_folder_full_path_ + "ker_" + suffix);
                        gt_exporter_.export_array_complex(A, debug_folder_full_path_ + "A_" + suffix);
                        gt_exporter_.export_array_complex(B, debug_folder_full_path_ + "B_" + suffix);
                    }

                    Gadgetron::grappa2d_convert_to_convolution_kernel(ker, kRO, kE1, oE1, convKer);
                    Gadgetron::grappa2d_image_domain_kernel(convKer, RO, E1, kIm);

                    if (!debug_folder_full_path_.empty())
                    {
                        gt_exporter_.export_array_complex(convKer, debug_folder_full_path_ + "convKer_" + suffix);
                        gt_exporter_.export_array_complex(kIm, debug_folder_full_path_ + "kIm_" + suffix);
                    }

                    hoNDArray<T> coilMap(RO, E1, dstCHA,
                        &(recon_obj.coil_map_(0, 0, 0, 0, n, s, slc)));
                    hoNDArray<T> unmixC(RO, E1, srcCHA,
                        &(recon_obj.unmixing_coeff_(0, 0, 0, 0, n, s, slc)));
                    hoNDArray<float> gFactor;

                    Gadgetron::grappa2d_unmixing_coeff(kIm, coilMap, (size_t)acceFactorE1_[e], unmixC, gFactor);
                    memcpy(&(recon_obj.gfactor_(0, 0, 0, 0, n, s, slc)), gFactor.begin(),
                        gFactor.get_number_of_bytes());

                    if (!debug_folder_full_path_.empty())
                    {
                        gt_exporter_.export_array_complex(coilMap, debug_folder_full_path_ + "coilMap_" + suffix);
                        gt_exporter_.export_array_complex(unmixC, debug_folder_full_path_ + "unmixC_" + suffix);
                        gt_exporter_.export_array(gFactor, debug_folder_full_path_ + "gFactor_" + suffix);
                    }

                    // ----------------------------------------
                    // grappa ai
                    // ----------------------------------------
                    models_[e][ii] = create_grappa_ai_model(ker);
                    bp::incref(models_[e][ii].ptr());
                }

                // -----------------------------------
            }
        }
    }

    void GenericReconCartesianGrappaAIGadget::perform_unwrapping(IsmrmrdReconBit &recon_bit, ReconObjType &recon_obj, size_t e)
    {
        size_t RO = recon_bit.data_.data_.get_size(0);
        size_t E1 = recon_bit.data_.data_.get_size(1);
        size_t E2 = recon_bit.data_.data_.get_size(2);
        size_t CHA = recon_bit.data_.data_.get_size(3);
        size_t N = recon_bit.data_.data_.get_size(4);
        size_t S = recon_bit.data_.data_.get_size(5);
        size_t SLC = recon_bit.data_.data_.get_size(6);

        hoNDArray<T> &src = recon_obj.ref_calib_;
        hoNDArray<T> &dst = recon_obj.ref_calib_dst_;

        size_t ref_RO = src.get_size(0);
        size_t ref_E1 = src.get_size(1);
        size_t ref_E2 = src.get_size(2);
        size_t srcCHA = src.get_size(3);
        size_t ref_N = src.get_size(4);
        size_t ref_S = src.get_size(5);
        size_t ref_SLC = src.get_size(6);

        size_t dstCHA = dst.get_size(3);

        size_t unmixingCoeff_CHA = recon_obj.unmixing_coeff_.get_size(3);

        size_t convkRO = recon_obj.kernel_.get_size(0);
        size_t convkE1 = recon_obj.kernel_.get_size(1);
        size_t convkE2 = recon_obj.kernel_.get_size(2);

        recon_obj.recon_res_.data_.create(RO, E1, E2, 1, N, S, SLC);
        recon_res_grappa_ai_[e].data_.create(RO, E1, E2, 1, N, S, SLC);

        /*if (!debug_folder_full_path_.empty())
        {
            std::stringstream os;
            os << "encoding_" << e;
            std::string suffix = os.str();
            gt_exporter_.export_array_complex(recon_bit.data_.data_, debug_folder_full_path_ + "data_src_" + suffix);
        }*/

        // compute aliased images
        data_recon_buf_.create(RO, E1, E2, CHA, N, S, SLC);

        if (E2 > 1)
        {
            Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_bit.data_.data_, complex_im_recon_buf_,
                data_recon_buf_);
        }
        else
        {
            Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_bit.data_.data_, complex_im_recon_buf_,
                data_recon_buf_);
        }

        // SNR unit scaling
        float effective_acce_factor(1), snr_scaling_ratio(1);
        this->compute_snr_scaling_factor(recon_bit, effective_acce_factor, snr_scaling_ratio);
        if (effective_acce_factor > 1)
        {
            // since the grappa in gadgetron is doing signal preserving scaling, to perserve noise level, we need this compensation factor
            double grappaKernelCompensationFactor = 1.0 / (acceFactorE1_[e] * acceFactorE2_[e]);
            Gadgetron::scal((float)(grappaKernelCompensationFactor * snr_scaling_ratio), complex_im_recon_buf_);

            if (this->verbose.value()) GDEBUG_STREAM(
                "GenericReconCartesianGrappaAIGadget, grappaKernelCompensationFactor*snr_scaling_ratio : "
                << grappaKernelCompensationFactor * snr_scaling_ratio);
        }

        if (!debug_folder_full_path_.empty())
        {
            std::stringstream os;
            os << "encoding_" << e;
            std::string suffix = os.str();
            gt_exporter_.export_array_complex(complex_im_recon_buf_, debug_folder_full_path_ + "aliasedIm_" + suffix);
        }

        // ------------------------------------------------
        // unwrapping with grappa
        // ------------------------------------------------

        long long num = N * S * SLC;
        long long ii;

        if (E2 > 1)
        {
            // to be implemented, 3D grappa ai
            recon_res_grappa_ai_[e].data_ = recon_obj.recon_res_.data_;
        }
        else
        {
            PythonFunction< hoNDArray<T> > apply_grappa_ai("grappa_ai", "apply_grappa_ai_model");

            size_t kRO = grappa_kSize_RO.value();
            size_t kNE1 = grappa_kSize_E1.value();

            std::vector<int> kE1, oE1;
            bool fitItself = false;
            size_t oNE1 = kernels_[e][0].get_size(4);
            if(oNE1 == acceFactorE1_[e]) fitItself = true;

            size_t convkRO, convkE1;
            Gadgetron::grappa2d_kerPattern(kE1, oE1, convkRO, convkE1, (size_t)acceFactorE1_[e], kRO, kNE1, fitItself);

            bool periodic_boundary_condition = true;

//#pragma omp parallel default(none) private(ii) shared(num, N, S, kRO, kE1, oE1, RO, E1, E2, srcCHA, convkRO, convkE1, convkE2, ref_N, ref_S, recon_obj, dstCHA, unmixingCoeff_CHA, e, recon_bit, periodic_boundary_condition, apply_grappa_ai) if(num>1)
            {
                hoNDArray<T> dataA;
                hoNDArray<unsigned short> dataAInd;
                hoNDArray<T> recon;
                hoNDArray<T> res_grappa, im_grappa, im_combined_grappa;
                hoNDArray<T> res_ai, im_ai, im_combined_ai;

//#pragma omp for
                for (ii = 0; ii < num; ii++)
                {
                    size_t slc = ii / (N * S);
                    size_t s = (ii - slc * N * S) / N;
                    size_t n = ii - slc * N * S - s * N;

                    // prepare recon
                    hoNDArray<T> data;
                    T *pData = &(recon_bit.data_.data_(0, 0, 0, 0, n, s, slc));
                    data.create(RO, E1, srcCHA, pData);

                    size_t usedN = n;
                    if (n >= ref_N) usedN = ref_N - 1;

                    size_t usedS = s;
                    if (s >= ref_S) usedS = ref_S - 1;

                    size_t ref_ii = usedN + usedS * ref_N + slc * ref_N*ref_S;

                    Gadgetron::grappa2d_prepare_recon(data, kRO, kE1, oE1, periodic_boundary_condition, dataA, dataAInd);

                    hoNDArray<T> coilMap(RO, E1, dstCHA, &(recon_obj.coil_map_(0, 0, 0, 0, usedN, usedS, slc)));

                    // grappa recon
                    res_grappa = data;
                    Gadgetron::grappa2d_perform_recon(dataA, kernels_[e][ref_ii], dataAInd, oE1, RO, E1, res_grappa);

                    Gadgetron::hoNDFFT<float>::instance()->ifft2c(res_grappa, im_grappa);
                    Gadgetron::coil_combine(im_grappa, coilMap, 2, im_combined_grappa);

                    memcpy(&(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc)), im_combined_grappa.begin(), im_combined_grappa.get_number_of_bytes());

                    // grappa ai recon
                    im_ai = im_grappa;
                    recon = apply_grappa_ai(dataA, models_[e][ref_ii]);
                    res_ai = data;
                    Gadgetron::grappa2d_fill_reconed_kspace(dataAInd, recon, oE1, RO, E1, res_ai);
                    Gadgetron::hoNDFFT<float>::instance()->ifft2c(res_ai, im_ai);
                    Gadgetron::coil_combine(im_ai, coilMap, 2, im_combined_ai);

                    memcpy(&(recon_res_grappa_ai_[e].data_(0, 0, 0, 0, n, s, slc)), im_combined_ai.begin(), im_combined_ai.get_number_of_bytes());
                }
            }
        }

        if (!debug_folder_full_path_.empty())
        {
            std::stringstream os;
            os << "encoding_" << e;
            std::string suffix = os.str();
            gt_exporter_.export_array_complex(recon_res_grappa_ai_[e].data_,
                debug_folder_full_path_ + "unwrappedIm_grappa_ai_" + suffix);
        }
    }

    int GenericReconCartesianGrappaAIGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "GenericReconCartesianGrappaAIGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (flags != 0)
        {
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GenericReconCartesianGrappaAIGadget)
}
