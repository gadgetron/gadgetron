
#include "GenericReconCartesianSpiritGadget.h"
#include "mri_core_spirit.h"
#include "hoNDArray_reductions.h"
#include "hoSPIRIT2DOperator.h"
#include "hoLsqrSolver.h"
#include "mri_core_grappa.h"

namespace Gadgetron {

    GenericReconCartesianSpiritGadget::GenericReconCartesianSpiritGadget() : BaseClass()
    {
    }

    GenericReconCartesianSpiritGadget::~GenericReconCartesianSpiritGadget()
    {
    }

    int GenericReconCartesianSpiritGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        // -------------------------------------------------

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        size_t NE = h.encoding.size();
        num_encoding_spaces_ = NE;
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        recon_obj_.resize(NE);

        // -------------------------------------------------
        // check the parameters
        if(this->spirit_iter_max.value()==0)
        {
            size_t acceFactor = acceFactorE1_[0] * acceFactorE2_[0];

            if(acceFactor>=6)
            {
                this->spirit_iter_max.value(150);
            }
            else if (acceFactor >= 5)
            {
                this->spirit_iter_max.value(120);
            }
            else if (acceFactor >= 4)
            {
                this->spirit_iter_max.value(100);
            }
            else if (acceFactor >= 3)
            {
                this->spirit_iter_max.value(60);
            }
            else
            {
                this->spirit_iter_max.value(50);
            }

            GDEBUG_STREAM("spirit_iter_max: " << this->spirit_iter_max.value());
        }

        if (this->spirit_iter_thres.value()<FLT_EPSILON)
        {
            this->spirit_iter_thres.value(0.0015);
            GDEBUG_STREAM("spirit_iter_thres: " << this->spirit_iter_thres.value());
        }

        if (this->spirit_reg_lamda.value()<FLT_EPSILON)
        {
            if(acceFactorE2_[0]>1)
            {
                this->spirit_reg_lamda.value(0.01);
            }
            else
            {
                this->spirit_reg_lamda.value(0.005);
            }
            GDEBUG_STREAM("spirit_reg_lamda: " << this->spirit_reg_lamda.value());
        }

        return GADGET_OK;
    }

    int GenericReconCartesianSpiritGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("GenericReconCartesianSpiritGadget::process"); }

        process_called_times_++;

        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
        }

        GadgetContainerMessage<std::vector<Core::Waveform>>* wav = AsContainerMessage<std::vector<Core::Waveform>>(m1->cont());
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
            os << "_encoding_" << e;

            GDEBUG_CONDITION_STREAM(verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(verbose.value(), "======================================================================");

            // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_, debug_folder_full_path_ + "data" + os.str()); }

            if (recon_bit_->rbit_[e].ref_)
            {
                // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_bit_->rbit_[e].ref_->data_, debug_folder_full_path_ + "ref" + os.str()); }

                // after this step, the recon_obj_[e].ref_calib_ and recon_obj_[e].ref_coil_map_ are set

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianSpiritGadget::make_ref_coil_map"); }
                this->make_ref_coil_map(*recon_bit_->rbit_[e].ref_,*recon_bit_->rbit_[e].data_.data_.get_dimensions(), recon_obj_[e].ref_calib_, recon_obj_[e].ref_coil_map_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // if (!debug_folder_full_path_.empty()) { this->gt_exporter_.export_array_complex(recon_obj_[e].ref_calib_, debug_folder_full_path_ + "ref_calib" + os.str()); }
                // if (!debug_folder_full_path_.empty()) { this->gt_exporter_.export_array_complex(recon_obj_[e].ref_coil_map_, debug_folder_full_path_ + "ref_coil_map" + os.str()); }

                // ----------------------------------------------------------

                // after this step, coil map is computed and stored in recon_obj_[e].coil_map_
                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianSpiritGadget::perform_coil_map_estimation"); }
                this->perform_coil_map_estimation(recon_obj_[e].ref_coil_map_, recon_obj_[e].coil_map_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // if (!debug_folder_full_path_.empty()) { this->gt_exporter_.export_array_complex(recon_obj_[e].coil_map_, debug_folder_full_path_ + "coil_map_" + os.str()); }

                // ---------------------------------------------------------------

                // after this step, recon_obj_[e].kernel_, recon_obj_[e].kernelIm_ or recon_obj_[e].kernelIm3D_ are filled
                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianSpiritGadget::perform_calib"); }
                this->perform_calib(recon_bit_->rbit_[e], recon_obj_[e], e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                // recon_bit_->rbit_[e].ref_ = Core::none;
            }

            if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0)
            {
                // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_, debug_folder_full_path_ + "data_before_unwrapping" + os.str()); }

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianSpiritGadget::perform_unwrapping"); }
                this->perform_unwrapping(recon_bit_->rbit_[e], recon_obj_[e], e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianSpiritGadget::compute_image_header"); }
                this->compute_image_header(recon_bit_->rbit_[e], recon_obj_[e].recon_res_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                if (wav) this->set_wave_form_to_image_array(*wav->getObjectPtr(), recon_obj_[e].recon_res_);
                recon_obj_[e].recon_res_.acq_headers_ = recon_bit_->rbit_[e].data_.headers_;

                // ---------------------------------------------------------------

                // if (!debug_folder_full_path_.empty()) { this->gt_exporter_.export_array_complex(recon_obj_[e].recon_res_.data_, debug_folder_full_path_ + "recon_res" + os.str()); }

                if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianSpiritGadget::send_out_image_array"); }
                this->send_out_image_array(recon_obj_[e].recon_res_, e, image_series.value() + ((int)e + 1), GADGETRON_IMAGE_REGULAR);
                if (perform_timing.value()) { gt_timer_.stop(); }
            }

//            recon_bit_->rbit_[e].ref_->clear();
            recon_bit_->rbit_[e].ref_ = Core::none;
            recon_obj_[e].recon_res_.data_.clear();
            recon_obj_[e].recon_res_.headers_.clear();
            recon_obj_[e].recon_res_.meta_.clear();
        }

        m1->release();

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    void GenericReconCartesianSpiritGadget::perform_calib(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            size_t RO = recon_bit.data_.data_.get_size(0);
            size_t E1 = recon_bit.data_.data_.get_size(1);
            size_t E2 = recon_bit.data_.data_.get_size(2);

            hoNDArray< std::complex<float> >& src = recon_obj.ref_calib_;
            hoNDArray< std::complex<float> >& dst = recon_obj.ref_calib_;

            size_t ref_RO = src.get_size(0);
            size_t ref_E1 = src.get_size(1);
            size_t ref_E2 = src.get_size(2);
            size_t srcCHA = src.get_size(3);
            size_t ref_N = src.get_size(4);
            size_t ref_S = src.get_size(5);
            size_t ref_SLC = src.get_size(6);

            size_t dstCHA = dst.get_size(3);

            recon_obj.gfactor_.create(RO, E1, E2, 1, ref_N, ref_S, ref_SLC);

            if (acceFactorE1_[e] > 1 || acceFactorE2_[e] > 1)
            {
                // allocate buffer for kernels
                size_t kRO = spirit_kSize_RO.value();
                size_t kE1 = spirit_kSize_E1.value();
                size_t kE2 = spirit_kSize_E2.value();

                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit, kRO : " << kRO);
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit, kE1 : " << kE1);
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit, kE2 : " << kE2);

                size_t convKRO = 2 * kRO - 1;
                size_t convKE1 = 2 * kE1 - 1;
                size_t convKE2 = 2 * kE2 - 1;

                if (E2 > 1)
                {
                    recon_obj.kernel_.create(convKRO, convKE1, convKE2, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);
                    recon_obj.kernelIm3D_.create(convKE1, convKE2, srcCHA, dstCHA, RO, ref_N, ref_S, ref_SLC);
                    Gadgetron::clear(recon_obj.kernelIm3D_);
                }
                else
                {
                    recon_obj.kernel_.create(convKRO, convKE1, 1, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);
                    recon_obj.kernelIm2D_.create(RO, E1, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);
                    Gadgetron::clear(recon_obj.kernelIm2D_);
                }

                Gadgetron::clear(recon_obj.kernel_);

                long long num = ref_N*ref_S*ref_SLC;

                double reg_lamda = this->spirit_reg_lamda.value();
                double over_determine_ratio = this->spirit_calib_over_determine_ratio.value();

                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit, reg_lamda : " << reg_lamda);
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit, over_determine_ratio : " << over_determine_ratio);

                long long ii;

#pragma omp parallel for default(none) private(ii) shared(src, dst, recon_obj, e, num, ref_N, ref_S, ref_RO, ref_E1, ref_E2, RO, E1, E2, dstCHA, srcCHA, convKRO, convKE1, convKE2, kRO, kE1, kE2, reg_lamda, over_determine_ratio) if(num>1)
                for (ii = 0; ii < num; ii++)
                {
                    size_t slc = ii / (ref_N*ref_S);
                    size_t s = (ii - slc*ref_N*ref_S) / (ref_N);
                    size_t n = ii - slc*ref_N*ref_S - s*ref_N;

                    std::stringstream os;
                    os << "n" << n << "_s" << s << "_slc" << slc << "_encoding_" << e;
                    std::string suffix = os.str();

                    std::complex<float>* pSrc = &(src(0, 0, 0, 0, n, s, slc));
                    hoNDArray< std::complex<float> > ref_src(ref_RO, ref_E1, ref_E2, srcCHA, pSrc);

                    std::complex<float>* pDst = &(dst(0, 0, 0, 0, n, s, slc));
                    hoNDArray< std::complex<float> > ref_dst(ref_RO, ref_E1, ref_E2, dstCHA, pDst);

                    // -----------------------------------

                    if (E2 > 1)
                    {
                        hoNDArray< std::complex<float> > convKer(convKRO, convKE1, convKE2, srcCHA, dstCHA, &(recon_obj.kernel_(0, 0, 0, 0, 0, n, s, slc)));
                        hoNDArray< std::complex<float> > kIm(convKE1, convKE2, srcCHA, dstCHA, RO, &(recon_obj.kernelIm3D_(0, 0, 0, 0, 0, n, s, slc)));

                        Gadgetron::spirit3d_calib_convolution_kernel(ref_src, ref_dst, reg_lamda, over_determine_ratio, kRO, kE1, kE2, 1, 1, 1, convKer, true);
                        Gadgetron::spirit3d_kspace_image_domain_kernel(convKer, RO, kIm);
                    }
                    else
                    {
                        hoNDArray< std::complex<float> > acsSrc(ref_RO, ref_E1, srcCHA, const_cast< std::complex<float>*>(ref_src.begin()));
                        hoNDArray< std::complex<float> > acsDst(ref_RO, ref_E1, dstCHA, const_cast< std::complex<float>*>(ref_dst.begin()));

                        hoNDArray< std::complex<float> > convKer(convKRO, convKE1, srcCHA, dstCHA, &(recon_obj.kernel_(0, 0, 0, 0, 0, n, s, slc)));
                        hoNDArray< std::complex<float> > kIm(RO, E1, srcCHA, dstCHA, &(recon_obj.kernelIm2D_(0, 0, 0, 0, n, s, slc)));

                        Gadgetron::spirit2d_calib_convolution_kernel(acsSrc, acsDst, reg_lamda, kRO, kE1, 1, 1, convKer, true);
                        Gadgetron::spirit2d_image_domain_kernel(convKer, RO, E1, kIm);

                        hoNDArray< std::complex<float> > convKerTmp, kImTmp;

                        Gadgetron::grappa2d_calib_convolution_kernel(acsSrc, acsDst, (size_t)acceFactorE1_[e], 0.005, 5, 4, convKerTmp);
                        Gadgetron::grappa2d_image_domain_kernel(convKerTmp, RO, E1, kImTmp);

                        hoNDArray< std::complex<float> > coilMap(RO, E1, dstCHA, &(recon_obj.coil_map_(0, 0, 0, 0, n, s, slc)));
                        hoNDArray< std::complex<float> > unmixC;
                        hoNDArray<float> gFactor;

                        Gadgetron::grappa2d_unmixing_coeff(kImTmp, coilMap, (size_t)acceFactorE1_[e], unmixC, gFactor);
                        memcpy(&(recon_obj.gfactor_(0, 0, 0, 0, n, s, slc)), gFactor.begin(), gFactor.get_number_of_bytes());
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianSpiritGadget::perform_calib(...) ... ");
        }
    }

    void GenericReconCartesianSpiritGadget::perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            typedef std::complex<float> T;

            size_t RO = recon_bit.data_.data_.get_size(0);
            size_t E1 = recon_bit.data_.data_.get_size(1);
            size_t E2 = recon_bit.data_.data_.get_size(2);
            size_t dstCHA = recon_bit.data_.data_.get_size(3);
            size_t N = recon_bit.data_.data_.get_size(4);
            size_t S = recon_bit.data_.data_.get_size(5);
            size_t SLC = recon_bit.data_.data_.get_size(6);

            hoNDArray< std::complex<float> >& src = recon_obj.ref_calib_;

            size_t ref_RO = src.get_size(0);
            size_t ref_E1 = src.get_size(1);
            size_t ref_E2 = src.get_size(2);
            size_t srcCHA = src.get_size(3);
            size_t ref_N = src.get_size(4);
            size_t ref_S = src.get_size(5);
            size_t ref_SLC = src.get_size(6);

            size_t convkRO = recon_obj.kernel_.get_size(0);
            size_t convkE1 = recon_obj.kernel_.get_size(1);
            size_t convkE2 = recon_obj.kernel_.get_size(2);

            recon_obj.recon_res_.data_.create(RO, E1, E2, 1, N, S, SLC);
            Gadgetron::clear(recon_obj.recon_res_.data_);
            recon_obj.full_kspace_ = recon_bit.data_.data_;
            Gadgetron::clear(recon_obj.full_kspace_);

            std::stringstream os;
            os << "encoding_" << e;
            std::string suffix = os.str();

            // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_bit.data_.data_, debug_folder_full_path_ + "data_src_" + suffix); }

            // ------------------------------------------------------------------
            // compute effective acceleration factor
            // ------------------------------------------------------------------
            float effective_acce_factor(1), snr_scaling_ratio(1);
            this->compute_snr_scaling_factor(recon_bit, effective_acce_factor, snr_scaling_ratio);
            if (effective_acce_factor > 1)
            {
                Gadgetron::scal(snr_scaling_ratio, recon_bit.data_.data_);
            }

            Gadgetron::GadgetronTimer timer(false);

            // ------------------------------------------------------------------
            // compute the reconstruction
            // ------------------------------------------------------------------
            long long num = N*S*SLC;
            long long ii;

            if(this->acceFactorE1_[e]<=1 && this->acceFactorE2_[e]<=1)
            {
                recon_obj.full_kspace_ = recon_bit.data_.data_;
            }
            else
            {
                hoNDArray< std::complex<float> >& kspace = recon_bit.data_.data_;
                hoNDArray< std::complex<float> >& res = recon_obj.full_kspace_;
                size_t iter_max = this->spirit_iter_max.value();
                double iter_thres = this->spirit_iter_thres.value();
                bool print_iter = this->spirit_print_iter.value();

                size_t RO_recon_size = 32; // every 32 images were computed together

                GDEBUG_CONDITION_STREAM(this->verbose.value(), "iter_max : " << iter_max);
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "iter_thres : " << iter_thres);
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "print_iter : " << print_iter);
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "RO_recon_size : " << RO_recon_size);

                if (E2 > 1)
                {
                    // 3D recon

                    hoNDArray<T> kspaceIfftRO(RO, E1, E2, srcCHA);
                    hoNDArray<T> kspaceIfftROPermuted(E1, E2, srcCHA, RO);
                    hoNDArray<T> kIm(E1, E2, srcCHA, dstCHA, RO_recon_size, 1, 1);
                    hoNDArray<T> res_ro_recon(E1, E2, 1, dstCHA, RO_recon_size, 1, 1);

                    for (ii = 0; ii < num; ii++)
                    {
                        size_t slc = ii / (N*S);
                        size_t s = (ii - slc*N*S) / N;
                        size_t n = ii - slc*N*S - s*N;

                        GDEBUG_CONDITION_STREAM(this->verbose.value(), "3D recon, [n s slc] : [" << n << " " << s << " " << slc << "]");

                        std::stringstream os;
                        os << "encoding_" << e << "_n" << n << "_s" << s << "_slc" << slc;
                        std::string suffix_3D = os.str();

                        size_t ro, e1, e2, scha, dcha;

                        // ------------------------------------------------------
                        // check whether the kspace is undersampled
                        // ------------------------------------------------------
                        bool undersampled = false;
                        for (e2= 0; e2< E2; e2++)
                        {
                            for (e1 = 0; e1 < E1; e1++)
                            {
                                if ((std::abs(kspace(RO / 2, e1, e2, srcCHA - 1, n, s, slc)) == 0)
                                    && (std::abs(kspace(RO / 2, e1, e2, 0, n, s, slc)) == 0))
                                {
                                    undersampled = true;
                                    break;
                                }
                            }
                        }

                        size_t kerN = (n<ref_N) ? n : ref_N-1;
                        size_t kerS = (s<ref_S) ? s : ref_S - 1;

                        // ---------------------------------------------------------------------
                        // permute the kspace
                        // ---------------------------------------------------------------------
                        std::complex<float>* pKSpace = &(kspace(0, 0, 0, 0, n, s, slc));
                        hoNDArray< std::complex<float> > kspace3D(RO, E1, E2, srcCHA, pKSpace);

                        if (this->perform_timing.value()) timer.start("SPIRIT linear 3D, ifft1c along RO ... ");
                        Gadgetron::hoNDFFT<float>::instance()->ifft1c(kspace3D, kspaceIfftRO);
                        if (this->perform_timing.value()) timer.stop();
                        // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(kspaceIfftRO, debug_folder_full_path_ + "kspaceIfftRO_" + suffix_3D); }

                        if (this->perform_timing.value()) timer.start("SPIRIT linear 3D, permute along RO ... ");
                        std::complex<float>* pKspaceRO = kspaceIfftRO.begin();
                        std::complex<float>* pKspacePermutedRO = kspaceIfftROPermuted.begin();
                        for (scha = 0; scha < srcCHA; scha++)
                        {
                            for (e2 = 0; e2 < E2; e2++)
                            {
                                for (e1 = 0; e1 < E1; e1++)
                                {
                                    size_t ind = e1*RO + e2*RO*E1 + scha*RO*E1*E2;
                                    size_t ind_dst = e1 + e2*E1 + scha*E1*E2;

                                    for (ro = 0; ro < RO; ro++)
                                    {
                                        pKspacePermutedRO[ind_dst + ro*E1*E2*srcCHA] = pKspaceRO[ro + ind];
                                    }
                                }
                            }
                        }
                        if (this->perform_timing.value()) timer.stop();
                        // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(kspaceIfftROPermuted, debug_folder_full_path_ + "kspaceIfftROPermuted_" + suffix_3D); }

                        // ---------------------------------------------------------------------
                        // get the kspace for recon
                        // ---------------------------------------------------------------------
                        hoNDArray< std::complex<float> > kspace3D_recon(E1, E2, 1, srcCHA, RO, kspaceIfftROPermuted.begin());

                        // ---------------------------------------------------------------------
                        // get spirit kernel for recon
                        // ---------------------------------------------------------------------
                        std::complex<float>* pKer = &(recon_obj.kernelIm3D_(0, 0, 0, 0, 0, kerN, kerS, slc));
                        hoNDArray< std::complex<float> > kIm3D_recon(convkE1, convkE2, srcCHA, dstCHA, RO, pKer);

                        // ---------------------------------------------------------------------
                        // get the array to store results
                        // ---------------------------------------------------------------------
                        std::complex<float>* pRes = &(res(0, 0, 0, 0, n, s, slc));
                        hoNDArray< std::complex<float> > res_recon(RO, E1, E2, dstCHA, pRes);

                        // ---------------------------------------------------------------------
                        // perform recon along RO
                        // ---------------------------------------------------------------------
                        size_t start_ro = 0;
                        for (start_ro = 0; start_ro < RO; start_ro += RO_recon_size)
                        {
                            size_t end_ro = start_ro + RO_recon_size - 1;
                            if (end_ro > RO) end_ro = RO - 1;
                            size_t num = end_ro - start_ro + 1;

                            GDEBUG_CONDITION_STREAM(this->verbose.value(), "3D recon, start_ro - end_ro : " << start_ro << " - " << end_ro);

                            hoNDArray< std::complex<float> > kspace3D_recon_ro(E1, E2, 1, srcCHA, num, 1, 1, kspace3D_recon.begin() + start_ro*E1*E2*srcCHA);
                            hoNDArray< std::complex<float> > kIm3D_recon_ro(convkE1, convkE2, srcCHA, dstCHA, num, 1, 1, kIm3D_recon.begin() + start_ro*convkE1*convkE2*srcCHA*dstCHA);

                            std::stringstream os_ro;
                            os_ro << "encoding_" << e << "_n" << n << "_s" << s << "_slc" << slc << "_ro" << start_ro << "_" << end_ro;
                            std::string suffix_3D_ro = os_ro.str();

                            // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(kspace3D_recon_ro, debug_folder_full_path_ + "kspace3D_recon_ro_" + suffix_3D_ro); }
                            // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(kIm3D_recon_ro, debug_folder_full_path_ + "kIm3D_recon_ro_" + suffix_3D_ro); }

                            if(num==RO_recon_size)
                            {
                                if (this->perform_timing.value()) timer.start("SPIRIT linear 3D, image domain kernel along E1 and E2 ... ");
                                Gadgetron::spirit3d_image_domain_kernel(kIm3D_recon_ro, E1, E2, kIm);
                                if (this->perform_timing.value()) timer.stop();
                                // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(kIm, debug_folder_full_path_ + "kIm_" + suffix_3D_ro); }

                                if (this->perform_timing.value()) timer.start("SPIRIT linear 3D, linear unwrapping ... ");
                                this->perform_spirit_unwrapping(kspace3D_recon_ro, kIm, res_ro_recon);
                                if (this->perform_timing.value()) timer.stop();
                            }
                            else
                            {
                                if (this->perform_timing.value()) timer.start("SPIRIT linear 3D, image domain kernel along E1 and E2, last ... ");
                                hoNDArray< std::complex<float> > kIm_last(E1, E2, srcCHA, dstCHA, num);
                                Gadgetron::spirit3d_image_domain_kernel(kIm3D_recon_ro, E1, E2, kIm_last);
                                if (this->perform_timing.value()) timer.stop();
                                // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(kIm_last, debug_folder_full_path_ + "kIm_last_" + suffix_3D_ro); }

                                if (this->perform_timing.value()) timer.start("SPIRIT linear 3D, linear unwrapping ... ");
                                this->perform_spirit_unwrapping(kspace3D_recon_ro, kIm_last, res_ro_recon);
                                if (this->perform_timing.value()) timer.stop();
                            }

                            // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(res_ro_recon, debug_folder_full_path_ + "res_ro_recon_" + suffix_3D_ro); }

                            // ---------------------------------------------------------------------
                            // copy the recon results
                            // ---------------------------------------------------------------------
                            std::complex<float>* pResRO = res_ro_recon.begin();
                            for (dcha = 0; dcha < dstCHA; dcha++)
                            {
                                for (e2 = 0; e2 < E2; e2++)
                                {
                                    for (e1 = 0; e1 < E1; e1++)
                                    {
                                        for (ro = 0; ro < num; ro++)
                                        {
                                            pRes[ro + start_ro + e1*RO + e2*RO*E1 + dcha*RO*E1*E2] = pResRO[e1 + e2*E1 + dcha*E1*E2 + ro*E1*E2*dstCHA];
                                        }
                                    }
                                }
                            }
                        }

                        // ---------------------------------------------------------------------
                        // go back to kspace for RO
                        // ---------------------------------------------------------------------
                        if (this->perform_timing.value()) timer.start("SPIRIT linear 3D, fft along RO for res ... ");
                        Gadgetron::hoNDFFT<float>::instance()->fft1c(res_recon);
                        if (this->perform_timing.value()) timer.stop();
                        // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(res_recon, debug_folder_full_path_ + "res_recon_" + suffix_3D); }
                    }
                }
                else
                {
                    if (this->perform_timing.value()) timer.start("SPIRIT 2D, linear unwrapping ... ");
                    this->perform_spirit_unwrapping(kspace, recon_obj.kernelIm2D_, res);
                    if (this->perform_timing.value()) timer.stop();

                    // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(res, debug_folder_full_path_ + "res_spirit_2D_" + suffix); }
                }
            }

            // ---------------------------------------------------------------------
            // compute coil combined images
            // ---------------------------------------------------------------------
            if (this->perform_timing.value()) timer.start("SPIRIT linear, coil combination ... ");
            this->perform_spirit_coil_combine(recon_obj);
            if (this->perform_timing.value()) timer.stop();

            // if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_obj.recon_res_.data_, debug_folder_full_path_ + "unwrappedIm_" + suffix); }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianSpiritGadget::perform_unwrapping(...) ... ");
        }
    }

    void GenericReconCartesianSpiritGadget::perform_spirit_coil_combine(ReconObjType& recon_obj)
    {
        try
        {
            size_t RO = recon_obj.full_kspace_.get_size(0);
            size_t E1 = recon_obj.full_kspace_.get_size(1);
            size_t E2 = recon_obj.full_kspace_.get_size(2);
            size_t dstCHA = recon_obj.full_kspace_.get_size(3);
            size_t N = recon_obj.full_kspace_.get_size(4);
            size_t S = recon_obj.full_kspace_.get_size(5);
            size_t SLC = recon_obj.full_kspace_.get_size(6);

            if (E2>1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_obj.full_kspace_, complex_im_recon_buf_);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_obj.full_kspace_, complex_im_recon_buf_);
            }

            size_t num = N*S*SLC;
            long long ii;

#pragma omp parallel default(none) private(ii) shared(num, N, S, recon_obj, RO, E1, E2, dstCHA) if(num>1)
            {
                hoNDArray< std::complex<float> > complexImBuf(RO, E1, E2, dstCHA);

#pragma omp for 
                for (ii = 0; ii < num; ii++)
                {
                    size_t slc = ii / (N*S);
                    size_t s = (ii - slc*N*S) / N;
                    size_t n = ii - slc*N*S - s*N;

                    size_t coilMapN = n;
                    if (coilMapN >= recon_obj.coil_map_.get_size(5)) coilMapN = recon_obj.coil_map_.get_size(5) - 1;

                    size_t coilMapS = s;
                    if (coilMapS >= recon_obj.coil_map_.get_size(6)) coilMapS = recon_obj.coil_map_.get_size(6) - 1;

                    hoNDArray< std::complex<float> > complexIm(RO, E1, E2, dstCHA, &(complex_im_recon_buf_(0, 0, 0, 0, n, s, slc)));
                    hoNDArray< std::complex<float> > coilMap(RO, E1, E2, dstCHA, &(recon_obj.coil_map_(0, 0, 0, 0, coilMapN, coilMapS, slc)));
                    hoNDArray< std::complex<float> > combined(RO, E1, E2, 1, &(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc)));

                    Gadgetron::multiplyConj(complexIm, coilMap, complexImBuf);
                    Gadgetron::sum_over_dimension(complexImBuf, combined, 3);
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianSpiritGadget::perform_spirit_coil_combine(...) ... ");
        }
    }

    void GenericReconCartesianSpiritGadget::perform_spirit_unwrapping(hoNDArray< std::complex<float> >& kspace, hoNDArray< std::complex<float> >& kerIm, hoNDArray< std::complex<float> >& res)
    {
        try
        {
            size_t iter_max = this->spirit_iter_max.value();
            double iter_thres = this->spirit_iter_thres.value();
            bool print_iter = this->spirit_print_iter.value();

            size_t RO = kspace.get_size(0);
            size_t E1 = kspace.get_size(1);
            size_t E2 = kspace.get_size(2);
            size_t CHA = kspace.get_size(3);
            size_t N = kspace.get_size(4);
            size_t S = kspace.get_size(5);
            size_t SLC = kspace.get_size(6);

            size_t ref_N = kerIm.get_size(4);
            size_t ref_S = kerIm.get_size(5);

            long long num = N*S*SLC;
            long long ii;

            hoNDArray< std::complex<float> > ker_Shifted(kerIm);
            Gadgetron::hoNDFFT<float>::instance()->ifftshift2D(kerIm, ker_Shifted);

            hoNDArray< std::complex<float> > kspace_Shifted;
            kspace_Shifted = kspace;
            Gadgetron::hoNDFFT<float>::instance()->ifftshift2D(kspace, kspace_Shifted);

#ifdef USE_OMP
            int numThreads = (int)num;
            if (numThreads > omp_get_num_procs()) numThreads = omp_get_num_procs();
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "numThreads : " << numThreads);
#endif // USE_OMP

            std::vector<size_t> dim(3, 1);
            dim[0] = RO;
            dim[1] = E1;
            dim[2] = CHA;

#pragma omp parallel default(none) private(ii) shared(num, N, S, RO, E1, CHA, dim, ref_N, ref_S, kspace, res, kspace_Shifted, ker_Shifted, iter_max, iter_thres, print_iter) num_threads(numThreads) if(num>1) 
            {
                boost::shared_ptr< hoSPIRIT2DOperator< std::complex<float> > > oper(new hoSPIRIT2DOperator< std::complex<float> >(&dim));
                hoSPIRIT2DOperator< std::complex<float> >& spirit = *oper;
                spirit.use_non_centered_fft_ = true;
                spirit.no_null_space_ = false;

                if (ref_N == 1 && ref_S == 1)
                {
                    boost::shared_ptr<hoNDArray< std::complex<float> > > ker(new hoNDArray< std::complex<float> >(RO, E1, CHA, CHA, ker_Shifted.begin()));
                    spirit.set_forward_kernel(*ker, false);
                }

                hoLsqrSolver< std::complex<float> > cgSolver;
                cgSolver.set_tc_tolerance((float)iter_thres);
                cgSolver.set_max_iterations(iter_max);
                cgSolver.set_output_mode(print_iter ? hoLsqrSolver< std::complex<float> >::OUTPUT_VERBOSE : hoLsqrSolver< std::complex<float> >::OUTPUT_SILENT);
                cgSolver.set_encoding_operator(oper);

                hoNDArray< std::complex<float> > b(RO, E1, CHA);
                hoNDArray< std::complex<float> > unwarppedKSpace(RO, E1, CHA);

#pragma omp for 
                for (ii = 0; ii < num; ii++)
                {
                    size_t slc = ii / (N*S);
                    size_t s = (ii - slc*N*S) / N;
                    size_t n = ii - slc*N*S - s*N;

                    // check whether the kspace is undersampled
                    bool undersampled = false;
                    for (size_t e1 = 0; e1 < E1; e1++)
                    {
                        if ((std::abs(kspace(RO / 2, e1, 0, CHA - 1, n, s, slc)) == 0)
                            && (std::abs(kspace(RO / 2, e1, 0, 0, n, s, slc)) == 0))
                        {
                            undersampled = true;
                            break;
                        }
                    }

                    std::complex<float>* pKpaceShifted = &(kspace_Shifted(0, 0, 0, 0, n, s, slc));
                    std::complex<float>* pRes = &(res(0, 0, 0, 0, n, s, slc));

                    if (!undersampled)
                    {
                        memcpy(pRes, pKpaceShifted, sizeof(std::complex<float>)*RO*E1*CHA);
                        continue;
                    }

                    long long kernelN = n;
                    if (kernelN >= (long long)ref_N) kernelN = (long long)ref_N - 1;

                    long long kernelS = s;
                    if (kernelS >= (long long)ref_S) kernelS = (long long)ref_S - 1;

                    boost::shared_ptr< hoNDArray< std::complex<float> > > acq(new hoNDArray< std::complex<float> >(RO, E1, CHA, pKpaceShifted));
                    spirit.set_acquired_points(*acq);
                    cgSolver.set_x0(acq);

                    if (ref_N == 1 && ref_S == 1)
                    {
                        spirit.compute_righ_hand_side(*acq, b);
                        cgSolver.solve(&unwarppedKSpace, &b);
                    }
                    else
                    {
                        std::complex<float>* pKer = &(ker_Shifted(0, 0, 0, 0, kernelN, kernelS, slc));
                        boost::shared_ptr<hoNDArray< std::complex<float> > > ker(new hoNDArray< std::complex<float> >(RO, E1, CHA, CHA, pKer));
                        spirit.set_forward_kernel(*ker, false);

                        spirit.compute_righ_hand_side(*acq, b);
                        cgSolver.solve(&unwarppedKSpace, &b);
                    }

                    // restore the acquired points
                    spirit.restore_acquired_kspace(*acq, unwarppedKSpace);
                    memcpy(pRes, unwarppedKSpace.begin(), unwarppedKSpace.get_number_of_bytes());
                }
            }

            Gadgetron::hoNDFFT<float>::instance()->fftshift2D(res, kspace_Shifted);
            res = kspace_Shifted;
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianSpiritGadget::perform_spirit_unwrapping(...) ... ");
        }
    }

    GADGET_FACTORY_DECLARE(GenericReconCartesianSpiritGadget)
}
