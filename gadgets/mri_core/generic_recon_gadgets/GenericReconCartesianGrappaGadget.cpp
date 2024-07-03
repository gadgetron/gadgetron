
#include "GenericReconCartesianGrappaGadget.h"
#include "mri_core_grappa.h"
#include "hoNDArray_reductions.h"

/*
    The input is IsmrmrdReconData and output is single 2D or 3D ISMRMRD images

    If required, the gfactor map can be sent out

    If the  number of required destination channel is 1, the GrappaONE recon will be performed

    The image number computation logic is implemented in compute_image_number function, which can be overloaded
*/

namespace Gadgetron {

    GenericReconCartesianGrappaGadget::GenericReconCartesianGrappaGadget() : BaseClass() {
    }

    GenericReconCartesianGrappaGadget::~GenericReconCartesianGrappaGadget() {
    }

    int GenericReconCartesianGrappaGadget::process_config(ACE_Message_Block *mb) {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        // -------------------------------------------------

        ISMRMRD::IsmrmrdHeader h;
        try {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...) {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        size_t NE = h.encoding.size();
        num_encoding_spaces_ = NE;
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        recon_obj_.resize(NE);

        GDEBUG("PATHNAME %s 'n",this->context.paths.gadgetron_home.c_str());

        this->stream_ismrmrd_header(h);

        return GADGET_OK;
    }

    int GenericReconCartesianGrappaGadget::process(Gadgetron::GadgetContainerMessage<IsmrmrdReconData> *m1) {
        if (perform_timing.value()) { gt_timer_local_.start("GenericReconCartesianGrappaGadget::process"); }

        process_called_times_++;

        IsmrmrdReconData *recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_) {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size()
                                                                                            << " instead of "
                                                                                            << num_encoding_spaces_);
        }

        GadgetContainerMessage<std::vector<Core::Waveform>>* wav = AsContainerMessage<std::vector<Core::Waveform>>(m1->cont());
        if (wav)
        {
            if (verbose.value())
            {
                GDEBUG_STREAM("Incoming recon_bit with " << wav->getObjectPtr()->size() << " wave form samples ");
            }
        }

        if (verbose.value())
        {
            for(auto key : this->buffer_names_)
            {
                GDEBUG_STREAM("buffer_names_ has " << key.first << " - " << key.second);
            }
        }

        // for every encoding space
        for (size_t e = 0; e < recon_bit_->rbit_.size(); e++) {
            std::stringstream os;
            os << "_encoding_" << e << "_" << process_called_times_;

            GDEBUG_CONDITION_STREAM(verbose.value(),
                                    "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(verbose.value(),
                                    "======================================================================");

            // ---------------------------------------------------------------
            // export incoming data
            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_, debug_folder_full_path_ + "data" + os.str());
            }

            if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].data_.trajectory_)
            {
                if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0)
                {
                    gt_exporter_.export_array(*(recon_bit_->rbit_[e].data_.trajectory_), debug_folder_full_path_ + "data_traj" + os.str());
                }
            }

            // ---------------------------------------------------------------

            if (recon_bit_->rbit_[e].ref_) {
                this->stream_to_array_buffer(GENERIC_RECON_REF_KSPACE, recon_bit_->rbit_[e].ref_->data_);

                if (!debug_folder_full_path_.empty()) {
                    gt_exporter_.export_array_complex(recon_bit_->rbit_[e].ref_->data_,
                                                      debug_folder_full_path_ + "ref" + os.str());
                }

                if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].ref_->trajectory_) {
                    if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0) {
                        gt_exporter_.export_array(*(recon_bit_->rbit_[e].ref_->trajectory_),
                                                  debug_folder_full_path_ + "ref_traj" + os.str());
                    }
                }

                // ---------------------------------------------------------------

                // after this step, the recon_obj_[e].ref_calib_ and recon_obj_[e].ref_coil_map_ are set

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianGrappaGadget::make_ref_coil_map"); }
                this->make_ref_coil_map(*recon_bit_->rbit_[e].ref_, recon_bit_->rbit_[e].data_.data_.get_dimensions(),
                                        recon_obj_[e].ref_calib_, recon_obj_[e].ref_coil_map_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ----------------------------------------------------------
                // export prepared ref for calibration and coil map
                if (!debug_folder_full_path_.empty()) {
                    this->gt_exporter_.export_array_complex(recon_obj_[e].ref_calib_,
                                                            debug_folder_full_path_ + "ref_calib" + os.str());
                }

                if (!debug_folder_full_path_.empty()) {
                    this->gt_exporter_.export_array_complex(recon_obj_[e].ref_coil_map_,
                                                            debug_folder_full_path_ + "ref_coil_map" + os.str());
                }

                // ---------------------------------------------------------------
                // after this step, the recon_obj_[e].ref_calib_dst_ and recon_obj_[e].ref_coil_map_ are modified
                if (perform_timing.value()) {
                    gt_timer_.start("GenericReconCartesianGrappaGadget::prepare_down_stream_coil_compression_ref_data");
                }
                this->prepare_down_stream_coil_compression_ref_data(recon_obj_[e].ref_calib_,
                                                                    recon_obj_[e].ref_coil_map_,
                                                                    recon_obj_[e].ref_calib_dst_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                this->stream_to_array_buffer(GENERIC_RECON_REF_KSPACE_FOR_COILMAP, recon_obj_[e].ref_coil_map_);

                if (!debug_folder_full_path_.empty()) {
                    this->gt_exporter_.export_array_complex(recon_obj_[e].ref_calib_dst_,
                        debug_folder_full_path_ + "ref_calib_dst" + os.str());
                }

                if (!debug_folder_full_path_.empty()) {
                    this->gt_exporter_.export_array_complex(recon_obj_[e].ref_coil_map_,
                        debug_folder_full_path_ + "ref_coil_map_dst" + os.str());
                }

                // ---------------------------------------------------------------

                // after this step, coil map is computed and stored in recon_obj_[e].coil_map_
                if (perform_timing.value()) {
                    gt_timer_.start("GenericReconCartesianGrappaGadget::perform_coil_map_estimation");
                }
                this->perform_coil_map_estimation(recon_obj_[e].ref_coil_map_, recon_obj_[e].coil_map_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                // after this step, recon_obj_[e].kernel_, recon_obj_[e].kernelIm_, recon_obj_[e].unmixing_coeff_ are filled
                // gfactor is computed too
                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianGrappaGadget::perform_calib"); }
                this->perform_calib(recon_bit_->rbit_[e], recon_obj_[e], e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                recon_bit_->rbit_[e].ref_ = Core::none;
            }

            if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0) {

                this->stream_to_array_buffer(GENERIC_RECON_UNDERSAMPLED_KSPACE, recon_bit_->rbit_[e].data_.data_);

                if (!debug_folder_full_path_.empty()) {
                    gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_,
                                                      debug_folder_full_path_ + "data_before_unwrapping" + os.str());
                }

                if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].data_.trajectory_) {
                    if (recon_bit_->rbit_[e].data_.trajectory_->get_number_of_elements() > 0) {
                        gt_exporter_.export_array(*(recon_bit_->rbit_[e].data_.trajectory_),
                                                  debug_folder_full_path_ + "data_before_unwrapping_traj" + os.str());
                    }
                }

                // ---------------------------------------------------------------

                if (perform_timing.value()) {
                    gt_timer_.start("GenericReconCartesianGrappaGadget::perform_unwrapping");
                }
                this->perform_unwrapping(recon_bit_->rbit_[e], recon_obj_[e], e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                if (perform_timing.value()) {
                    gt_timer_.start("GenericReconCartesianGrappaGadget::compute_image_header");
                }
                this->compute_image_header(recon_bit_->rbit_[e], recon_obj_[e].recon_res_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------
                // pass down waveform
                if (wav) this->set_wave_form_to_image_array(*wav->getObjectPtr(), recon_obj_[e].recon_res_);
                recon_obj_[e].recon_res_.acq_headers_ = recon_bit_->rbit_[e].data_.headers_;

                // ---------------------------------------------------------------
                if (send_out_gfactor.value() && recon_obj_[e].gfactor_.get_number_of_elements() > 0 &&
                    (acceFactorE1_[e] * acceFactorE2_[e] > 1)) {
                    IsmrmrdImageArray res;
                    Gadgetron::real_to_complex(recon_obj_[e].gfactor_, res.data_);
                    res.headers_ = recon_obj_[e].recon_res_.headers_;
                    res.meta_ = recon_obj_[e].recon_res_.meta_;

                    if (!debug_folder_full_path_.empty()) {
                        gt_exporter_.export_array_complex(res.data_,
                                                        debug_folder_full_path_ + "gfactor_" + os.str());
                    }

                    if (perform_timing.value()) {
                        gt_timer_.start("GenericReconCartesianGrappaGadget::send_out_image_array, gfactor");
                    }
                    this->send_out_image_array(res, e, image_series.value() + 10 * ((int) e + 2),
                                               GADGETRON_IMAGE_GFACTOR);
                    if (perform_timing.value()) { gt_timer_.stop(); }
                }

                // ---------------------------------------------------------------
                if (send_out_snr_map.value()) {
                    hoNDArray<std::complex<float> > snr_map;

                    if (calib_mode_[e] == Gadgetron::ISMRMRD_noacceleration) {
                        snr_map = recon_obj_[e].recon_res_.data_;
                    } else {
                        if (recon_obj_[e].gfactor_.get_number_of_elements() > 0) {
                            if (perform_timing.value()) { gt_timer_.start("compute SNR map array"); }
                            this->compute_snr_map(recon_obj_[e], snr_map);
                            if (perform_timing.value()) { gt_timer_.stop(); }
                        }
                    }

                    if (snr_map.get_number_of_elements() > 0) {
                        if (!debug_folder_full_path_.empty()) {
                            this->gt_exporter_.export_array_complex(snr_map,
                                                                    debug_folder_full_path_ + "snr_map" + os.str());
                        }

                        if (perform_timing.value()) { gt_timer_.start("send out gfactor array, snr map"); }

                        IsmrmrdImageArray res;
                        res.data_ = snr_map;
                        res.headers_ = recon_obj_[e].recon_res_.headers_;
                        res.meta_ = recon_obj_[e].recon_res_.meta_;
                        res.acq_headers_ = recon_bit_->rbit_[e].data_.headers_;

                        this->send_out_image_array(res, e,
                                                   image_series.value() + 100 * ((int) e + 3), GADGETRON_IMAGE_SNR_MAP);

                        if (perform_timing.value()) { gt_timer_.stop(); }
                    }
                }

                // ---------------------------------------------------------------

                if (!debug_folder_full_path_.empty()) {
                    this->gt_exporter_.export_array_complex(recon_obj_[e].recon_res_.data_,
                        debug_folder_full_path_ + "recon_res" + os.str());
                }

                this->stream_to_ismrmrd_image_buffer(GENERIC_RECON_COILMAP, recon_obj_[e].coil_map_, recon_obj_[e].recon_res_.headers_, recon_obj_[e].recon_res_.meta_);
                if (recon_obj_[e].gfactor_.get_number_of_elements() > 0) this->stream_to_ismrmrd_image_buffer(GENERIC_RECON_GFACTOR_MAP, recon_obj_[e].gfactor_, recon_obj_[e].recon_res_.headers_, recon_obj_[e].recon_res_.meta_);
                this->stream_to_ismrmrd_image_buffer(GENERIC_RECON_RECONED_COMPLEX_IMAGE, recon_obj_[e].recon_res_.data_, recon_obj_[e].recon_res_.headers_, recon_obj_[e].recon_res_.meta_);

                if (perform_timing.value()) {
                    gt_timer_.start("GenericReconCartesianGrappaGadget::send_out_image_array");
                }

                this->send_out_image_array(recon_obj_[e].recon_res_, e,
                    image_series.value() + ((int)e + 1), GADGETRON_IMAGE_REGULAR);
                if (perform_timing.value()) { gt_timer_.stop(); }
            }

            recon_obj_[e].recon_res_.data_.clear();
            recon_obj_[e].gfactor_.clear();
            recon_obj_[e].recon_res_.headers_.clear();
            recon_obj_[e].recon_res_.meta_.clear();
        }

        m1->release();

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    void GenericReconCartesianGrappaGadget::prepare_down_stream_coil_compression_ref_data(
            const hoNDArray<std::complex<float> > &ref_src, hoNDArray<std::complex<float> > &ref_coil_map,
            hoNDArray<std::complex<float> > &ref_dst, size_t e) {

        if (!downstream_coil_compression.value()) {
            GDEBUG_CONDITION_STREAM(verbose.value(), "Downstream coil compression is not prescribed ... ");
            ref_dst = ref_src;
            return;
        }

        if (downstream_coil_compression_thres.value() < 0 && downstream_coil_compression_num_modesKept.value() == 0) {
            GDEBUG_CONDITION_STREAM(verbose.value(),
                                    "Downstream coil compression is prescribed to use all input channels ... ");
            ref_dst = ref_src;
            return;
        }

        // determine how many channels to use
        size_t RO = ref_src.get_size(0);
        size_t E1 = ref_src.get_size(1);
        size_t E2 = ref_src.get_size(2);
        size_t CHA = ref_src.get_size(3);
        size_t N = ref_src.get_size(4);
        size_t S = ref_src.get_size(5);
        size_t SLC = ref_src.get_size(6);

        size_t recon_RO = ref_coil_map.get_size(0);
        size_t recon_E1 = ref_coil_map.get_size(1);
        size_t recon_E2 = ref_coil_map.get_size(2);

        std::complex<float> *pRef = const_cast< std::complex<float> * >(ref_src.begin());

        size_t dstCHA = CHA;
        if (downstream_coil_compression_num_modesKept.value() > 0 &&
            downstream_coil_compression_num_modesKept.value() <= CHA) {
            dstCHA = downstream_coil_compression_num_modesKept.value();
        } else {
            std::vector<float> E(CHA, 0);
            long long cha;

#pragma omp parallel default(none) private(cha) shared(RO, E1, E2, CHA, pRef, E)
            {
                hoNDArray<std::complex<float> > dataCha;
#pragma omp for
                for (cha = 0; cha < (long long) CHA; cha++) {
                    dataCha.create(RO, E1, E2, pRef + cha * RO * E1 * E2);
                    float v = Gadgetron::nrm2(dataCha);
                    E[cha] = v * v;
                }
            }

            for (cha = 1; cha < (long long) CHA; cha++) {
                if (std::abs(E[cha]) < downstream_coil_compression_thres.value() * std::abs(E[0])) {
                    break;
                }
            }

            dstCHA = cha;
        }

        GDEBUG_CONDITION_STREAM(verbose.value(),
                                "Downstream coil compression is prescribed to use " << dstCHA << " out of " << CHA
                                                                                    << " channels ...");

        if (dstCHA < CHA) {
            ref_dst.create(RO, E1, E2, dstCHA, N, S, SLC);
            hoNDArray<std::complex<float> > ref_coil_map_dst;
            ref_coil_map_dst.create(recon_RO, recon_E1, recon_E2, dstCHA, N, S, SLC);

            size_t slc, s, n;
            for (slc = 0; slc < SLC; slc++) {
                for (s = 0; s < S; s++) {
                    for (n = 0; n < N; n++) {
                        std::complex<float> *pDst = &(ref_dst(0, 0, 0, 0, n, s, slc));
                        const std::complex<float> *pSrc = &(ref_src(0, 0, 0, 0, n, s, slc));
                        memcpy(pDst, pSrc, sizeof(std::complex<float>) * RO * E1 * E2 * dstCHA);

                        pDst = &(ref_coil_map_dst(0, 0, 0, 0, n, s, slc));
                        pSrc = &(ref_coil_map(0, 0, 0, 0, n, s, slc));
                        memcpy(pDst, pSrc, sizeof(std::complex<float>) * recon_RO * recon_E1 * recon_E2 * dstCHA);
                    }
                }
            }

            ref_coil_map = ref_coil_map_dst;
        } else {
            ref_dst = ref_src;
        }

    }

    void
    GenericReconCartesianGrappaGadget::perform_calib(IsmrmrdReconBit &recon_bit, ReconObjType &recon_obj, size_t e) {

        size_t RO = recon_bit.data_.data_.get_size(0);
        size_t E1 = recon_bit.data_.data_.get_size(1);
        size_t E2 = recon_bit.data_.data_.get_size(2);

        hoNDArray<std::complex<float> > &src = recon_obj.ref_calib_;
        hoNDArray<std::complex<float> > &dst = recon_obj.ref_calib_dst_;

        size_t ref_RO = src.get_size(0);
        size_t ref_E1 = src.get_size(1);
        size_t ref_E2 = src.get_size(2);
        size_t srcCHA = src.get_size(3);
        size_t ref_N = src.get_size(4);
        size_t ref_S = src.get_size(5);
        size_t ref_SLC = src.get_size(6);

        size_t dstCHA = dst.get_size(3);

        recon_obj.unmixing_coeff_.create(RO, E1, E2, srcCHA, ref_N, ref_S, ref_SLC);
        recon_obj.gfactor_.create(RO, E1, E2, 1, ref_N, ref_S, ref_SLC);

        Gadgetron::clear(recon_obj.unmixing_coeff_);
        Gadgetron::clear(recon_obj.gfactor_);

        if (acceFactorE1_[e] <= 1 && acceFactorE2_[e] <= 1) {
            Gadgetron::conjugate(recon_obj.coil_map_, recon_obj.unmixing_coeff_);
        } else {
            // allocate buffer for kernels
            size_t kRO = grappa_kSize_RO.value();
            size_t kNE1 = grappa_kSize_E1.value();
            size_t kNE2 = grappa_kSize_E2.value();

            size_t convKRO(1), convKE1(1), convKE2(1);

            bool fitItself = this->downstream_coil_compression.value();

            if (E2 > 1) {
                std::vector<int> kE1, oE1;
                std::vector<int> kE2, oE2;
                grappa3d_kerPattern(kE1, oE1, kE2, oE2, convKRO, convKE1, convKE2, (size_t) acceFactorE1_[e],
                                    (size_t) acceFactorE2_[e], kRO, kNE1, kNE2, fitItself);
            } else {
                std::vector<int> kE1, oE1;
                Gadgetron::grappa2d_kerPattern(kE1, oE1, convKRO, convKE1, (size_t) acceFactorE1_[e], kRO, kNE1,
                                               fitItself);
                recon_obj.kernelIm_.create(RO, E1, 1, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);
            }

            recon_obj.kernel_.create(convKRO, convKE1, convKE2, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);

            Gadgetron::clear(recon_obj.kernel_);
            Gadgetron::clear(recon_obj.kernelIm_);

            long long num = ref_N * ref_S * ref_SLC;

            long long ii;

            // only allow this for loop openmp if num>1 and 2D recon
#pragma omp parallel for default(none) private(ii) shared(src, dst, recon_obj, e, num, ref_N, ref_S, ref_RO, ref_E1, ref_E2, RO, E1, E2, dstCHA, srcCHA, convKRO, convKE1, convKE2, kRO, kNE1, kNE2, fitItself) if(num>1)
            for (ii = 0; ii < num; ii++) {
                size_t slc = ii / (ref_N * ref_S);
                size_t s = (ii - slc * ref_N * ref_S) / (ref_N);
                size_t n = ii - slc * ref_N * ref_S - s * ref_N;

                std::stringstream os;
                os << "n" << n << "_s" << s << "_slc" << slc << "_encoding_" << e;
                std::string suffix = os.str();

                std::complex<float> *pSrc = &(src(0, 0, 0, 0, n, s, slc));
                hoNDArray<std::complex<float> > ref_src(ref_RO, ref_E1, ref_E2, srcCHA, pSrc);

                std::complex<float> *pDst = &(dst(0, 0, 0, 0, n, s, slc));
                hoNDArray<std::complex<float> > ref_dst(ref_RO, ref_E1, ref_E2, dstCHA, pDst);

                // -----------------------------------

                if (E2 > 1) {
                    hoNDArray<std::complex<float> > ker(convKRO, convKE1, convKE2, srcCHA, dstCHA,
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

                    //if (!debug_folder_full_path_.empty())
                    //{
                    //    gt_exporter_.export_array_complex(ker, debug_folder_full_path_ + "convKer3D_" + suffix);
                    //}

                    hoNDArray<std::complex<float> > coilMap(RO, E1, E2, dstCHA,
                                                            &(recon_obj.coil_map_(0, 0, 0, 0, n, s, slc)));
                    hoNDArray<std::complex<float> > unmixC(RO, E1, E2, srcCHA,
                                                           &(recon_obj.unmixing_coeff_(0, 0, 0, 0, n, s, slc)));
                    hoNDArray<float> gFactor(RO, E1, E2, 1, &(recon_obj.gfactor_(0, 0, 0, 0, n, s, slc)));
                    Gadgetron::grappa3d_unmixing_coeff(ker, coilMap, (size_t) acceFactorE1_[e],
                                                       (size_t) acceFactorE2_[e], unmixC, gFactor);

                    //if (!debug_folder_full_path_.empty())
                    //{
                    //    gt_exporter_.export_array_complex(unmixC, debug_folder_full_path_ + "unmixC_3D_" + suffix);
                    //}

                    //if (!debug_folder_full_path_.empty())
                    //{
                    //    gt_exporter_.export_array(gFactor, debug_folder_full_path_ + "gFactor_3D_" + suffix);
                    //}
                } else {
                    hoNDArray<std::complex<float> > acsSrc(ref_RO, ref_E1, srcCHA,
                                                           const_cast< std::complex<float> *>(ref_src.begin()));
                    hoNDArray<std::complex<float> > acsDst(ref_RO, ref_E1, dstCHA,
                                                           const_cast< std::complex<float> *>(ref_dst.begin()));

                    hoNDArray<std::complex<float> > convKer(convKRO, convKE1, srcCHA, dstCHA,
                                                            &(recon_obj.kernel_(0, 0, 0, 0, 0, n, s, slc)));
                    hoNDArray<std::complex<float> > kIm(RO, E1, srcCHA, dstCHA,
                                                        &(recon_obj.kernelIm_(0, 0, 0, 0, 0, n, s, slc)));

                    if (fitItself)
                    {
                        Gadgetron::grappa2d_calib_convolution_kernel(acsSrc, acsDst, (size_t)acceFactorE1_[e],
                            grappa_reg_lamda.value(), kRO, kNE1, convKer);
                    }
                    else
                    {
                        Gadgetron::grappa2d_calib_convolution_kernel(acsSrc, acsSrc, (size_t)acceFactorE1_[e],
                            grappa_reg_lamda.value(), kRO, kNE1, convKer);
                    }
                    Gadgetron::grappa2d_image_domain_kernel(convKer, RO, E1, kIm);

                    /*if (!debug_folder_full_path_.empty())
                    {
                        gt_exporter_.export_array_complex(convKer, debug_folder_full_path_ + "convKer_" + suffix);
                    }

                    if (!debug_folder_full_path_.empty())
                    {
                        gt_exporter_.export_array_complex(kIm, debug_folder_full_path_ + "kIm_" + suffix);
                    }*/

                    hoNDArray<std::complex<float> > coilMap(RO, E1, dstCHA,
                                                            &(recon_obj.coil_map_(0, 0, 0, 0, n, s, slc)));
                    hoNDArray<std::complex<float> > unmixC(RO, E1, srcCHA,
                                                           &(recon_obj.unmixing_coeff_(0, 0, 0, 0, n, s, slc)));
                    hoNDArray<float> gFactor;

                    Gadgetron::grappa2d_unmixing_coeff(kIm, coilMap, (size_t) acceFactorE1_[e], unmixC, gFactor);
                    memcpy(&(recon_obj.gfactor_(0, 0, 0, 0, n, s, slc)), gFactor.begin(),
                           gFactor.get_number_of_bytes());

                    // if (!debug_folder_full_path_.empty())
                    // {
                    //     gt_exporter_.export_array_complex(unmixC, debug_folder_full_path_ + "unmixC_" + suffix);
                    // }

                    // if (!debug_folder_full_path_.empty())
                    // {
                    //     gt_exporter_.export_array(gFactor, debug_folder_full_path_ + "gFactor_" + suffix);
                    // }
                }

                // -----------------------------------
            }
        }

    }

    void GenericReconCartesianGrappaGadget::perform_unwrapping(IsmrmrdReconBit &recon_bit, ReconObjType &recon_obj,
                                                               size_t e) {

        typedef std::complex<float> T;

        typedef std::complex<float> T;

        size_t RO = recon_bit.data_.data_.get_size(0);
        size_t E1 = recon_bit.data_.data_.get_size(1);
        size_t E2 = recon_bit.data_.data_.get_size(2);
        size_t dstCHA = recon_bit.data_.data_.get_size(3);
        size_t N = recon_bit.data_.data_.get_size(4);
        size_t S = recon_bit.data_.data_.get_size(5);
        size_t SLC = recon_bit.data_.data_.get_size(6);

        hoNDArray<std::complex<float> > &src = recon_obj.ref_calib_;

        size_t ref_RO = src.get_size(0);
        size_t ref_E1 = src.get_size(1);
        size_t ref_E2 = src.get_size(2);
        size_t srcCHA = src.get_size(3);
        size_t ref_N = src.get_size(4);
        size_t ref_S = src.get_size(5);
        size_t ref_SLC = src.get_size(6);

        size_t unmixingCoeff_CHA = recon_obj.unmixing_coeff_.get_size(3);

        size_t convkRO = recon_obj.kernel_.get_size(0);
        size_t convkE1 = recon_obj.kernel_.get_size(1);
        size_t convkE2 = recon_obj.kernel_.get_size(2);

        recon_obj.recon_res_.data_.create(RO, E1, E2, 1, N, S, SLC);

        if (!debug_folder_full_path_.empty()) {
            std::stringstream os;
            os << "encoding_" << e;
            std::string suffix = os.str();
            gt_exporter_.export_array_complex(recon_bit.data_.data_, debug_folder_full_path_ + "data_src_" + suffix);
        }

        // compute aliased images
        data_recon_buf_.create(RO, E1, E2, dstCHA, N, S, SLC);

        if (E2 > 1) {
            Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_bit.data_.data_, complex_im_recon_buf_,
                                                          data_recon_buf_);
        } else {
            Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_bit.data_.data_, complex_im_recon_buf_,
                                                          data_recon_buf_);
        }

        // SNR unit scaling
        float effective_acce_factor(1), snr_scaling_ratio(1);
        this->compute_snr_scaling_factor(recon_bit, effective_acce_factor, snr_scaling_ratio);
        if (effective_acce_factor > 1) {
            // since the grappa in gadgetron is doing signal preserving scaling, to preserve noise level, we need this compensation factor
            double grappaKernelCompensationFactor = 1.0 / (acceFactorE1_[e] * acceFactorE2_[e]);
            Gadgetron::scal((float) (grappaKernelCompensationFactor * snr_scaling_ratio), complex_im_recon_buf_);

            if (this->verbose.value()) GDEBUG_STREAM(
                    "GenericReconCartesianGrappaGadget, grappaKernelCompensationFactor*snr_scaling_ratio : "
                            << grappaKernelCompensationFactor * snr_scaling_ratio);
        }

        if (!debug_folder_full_path_.empty()) {
            std::stringstream os;
            os << "encoding_" << e;
            std::string suffix = os.str();
            gt_exporter_.export_array_complex(complex_im_recon_buf_, debug_folder_full_path_ + "aliasedIm_" + suffix);
        }

        // unwrapping

        long long num = N * S * SLC;

        long long ii;

#pragma omp parallel default(none) private(ii) shared(num, N, S, RO, E1, E2, srcCHA, convkRO, convkE1, convkE2, ref_N, ref_S, recon_obj, dstCHA, unmixingCoeff_CHA, e) if(num>1)
        {
#pragma omp for
            for (ii = 0; ii < num; ii++) {
                size_t slc = ii / (N * S);
                size_t s = (ii - slc * N * S) / N;
                size_t n = ii - slc * N * S - s * N;

                // combined channels
                T *pIm = &(complex_im_recon_buf_(0, 0, 0, 0, n, s, slc));

                size_t usedN = n;
                if (n >= ref_N) usedN = ref_N - 1;

                size_t usedS = s;
                if (s >= ref_S) usedS = ref_S - 1;

                T *pUnmix = &(recon_obj.unmixing_coeff_(0, 0, 0, 0, usedN, usedS, slc));

                T *pRes = &(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc));
                hoNDArray<std::complex<float> > res(RO, E1, E2, 1, pRes);

                hoNDArray<std::complex<float> > unmixing(RO, E1, E2, unmixingCoeff_CHA, pUnmix);
                hoNDArray<std::complex<float> > aliasedIm(RO, E1, E2,
                                                          ((unmixingCoeff_CHA <= srcCHA) ? unmixingCoeff_CHA : srcCHA),
                                                          1, pIm);
                Gadgetron::apply_unmix_coeff_aliased_image_3D(aliasedIm, unmixing, res);
            }
        }

        if (!debug_folder_full_path_.empty()) {
            std::stringstream os;
            os << "encoding_" << e;
            std::string suffix = os.str();
            gt_exporter_.export_array_complex(recon_obj.recon_res_.data_,
                                              debug_folder_full_path_ + "unwrappedIm_" + suffix);
        }

    }

    void GenericReconCartesianGrappaGadget::compute_snr_map(ReconObjType &recon_obj,
                                                            hoNDArray<std::complex<float> > &snr_map) {

        typedef std::complex<float> T;

        snr_map = recon_obj.recon_res_.data_;

        size_t RO = recon_obj.recon_res_.data_.get_size(0);
        size_t E1 = recon_obj.recon_res_.data_.get_size(1);
        size_t E2 = recon_obj.recon_res_.data_.get_size(2);
        size_t CHA = recon_obj.recon_res_.data_.get_size(3);
        size_t N = recon_obj.recon_res_.data_.get_size(4);
        size_t S = recon_obj.recon_res_.data_.get_size(5);
        size_t SLC = recon_obj.recon_res_.data_.get_size(6);

        size_t gN = recon_obj.gfactor_.get_size(4);
        size_t gS = recon_obj.gfactor_.get_size(5);

        GADGET_CHECK_THROW(recon_obj.gfactor_.get_size(0) == RO);
        GADGET_CHECK_THROW(recon_obj.gfactor_.get_size(1) == E1);
        GADGET_CHECK_THROW(recon_obj.gfactor_.get_size(2) == E2);
        GADGET_CHECK_THROW(recon_obj.gfactor_.get_size(3) == CHA);
        GADGET_CHECK_THROW(recon_obj.gfactor_.get_size(6) == SLC);

        size_t n, s, slc;
        for (slc = 0; slc < SLC; slc++) {
            for (s = 0; s < S; s++) {
                size_t usedS = s;
                if (usedS >= gS) usedS = gS - 1;

                for (n = 0; n < N; n++) {
                    size_t usedN = n;
                    if (usedN >= gN) usedN = gN - 1;

                    float *pG = &(recon_obj.gfactor_(0, 0, 0, 0, usedN, usedS, slc));
                    T *pIm = &(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc));
                    T *pSNR = &(snr_map(0, 0, 0, 0, n, s, slc));

                    for (size_t ii = 0; ii < RO * E1 * E2 * CHA; ii++) {
                        pSNR[ii] = pIm[ii] / pG[ii];
                    }
                }
            }
        }

    }



    GADGET_FACTORY_DECLARE(GenericReconCartesianGrappaGadget)
}
