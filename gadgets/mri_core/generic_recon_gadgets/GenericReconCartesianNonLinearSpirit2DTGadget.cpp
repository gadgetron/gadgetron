
#include "GenericReconCartesianNonLinearSpirit2DTGadget.h"
#include "hoSPIRIT2DTOperator.h"
#include "hoSPIRIT2DTDataFidelityOperator.h"
#include "hoWavelet2DTOperator.h"
#include "mri_core_spirit.h"
#include "mri_core_grappa.h"
#include "hoNDArray_reductions.h"
#include "hoGdSolver.h"
#include <boost/make_shared.hpp>

namespace Gadgetron {

    GenericReconCartesianNonLinearSpirit2DTGadget::GenericReconCartesianNonLinearSpirit2DTGadget() : BaseClass()
    {
    }

    GenericReconCartesianNonLinearSpirit2DTGadget::~GenericReconCartesianNonLinearSpirit2DTGadget()
    {
    }

    int GenericReconCartesianNonLinearSpirit2DTGadget::process_config(ACE_Message_Block* mb)
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

        // -------------------------------------------------
        // check the parameters
        if(this->spirit_nl_iter_max.value()==0)
        {
            this->spirit_nl_iter_max.value(15);
            GDEBUG_STREAM("spirit_iter_max: " << this->spirit_nl_iter_max.value());
        }

        if (this->spirit_nl_iter_thres.value()<FLT_EPSILON)
        {
            this->spirit_nl_iter_thres.value(0.004);
            GDEBUG_STREAM("spirit_nl_iter_thres: " << this->spirit_nl_iter_thres.value());
        }

        if (this->spirit_image_reg_lamda.value() < FLT_EPSILON)
        {
            if(this->spirit_reg_proximity_across_cha.value())
            {
                if(spirit_reg_estimate_noise_floor.value())
                {
                    this->spirit_image_reg_lamda.value(0.001);
                }
                else
                {
                    this->spirit_image_reg_lamda.value(0.0002);
                }
            }
            else
            {
                if(spirit_reg_estimate_noise_floor.value())
                {
                    this->spirit_image_reg_lamda.value(0.002);
                }
                else
                {
                    this->spirit_image_reg_lamda.value(0.00005);
                }
            }

            GDEBUG_STREAM("spirit_image_reg_lamda: " << this->spirit_image_reg_lamda.value());
        }

        if (this->spirit_reg_N_weighting_ratio.value() < FLT_EPSILON)
        {
            if(acceFactorE1_[0]<=5)
            {
                this->spirit_reg_N_weighting_ratio.value(10.0);
            }
            else
            {
                this->spirit_reg_N_weighting_ratio.value(20.0);
            }

            GDEBUG_STREAM("spirit_reg_N_weighting_ratio: " << this->spirit_reg_N_weighting_ratio.value());
        }

        return GADGET_OK;
    }

    void GenericReconCartesianNonLinearSpirit2DTGadget::perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
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

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_bit.data_.data_, debug_folder_full_path_ + "data_src_" + suffix); }

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
            if(this->acceFactorE1_[e]<=1 && this->acceFactorE2_[e]<=1)
            {
                recon_obj.full_kspace_ = recon_bit.data_.data_;
            }
            else
            {
                hoNDArray< std::complex<float> >& kspace = recon_bit.data_.data_;
                hoNDArray< std::complex<float> >& res = recon_obj.full_kspace_;
                hoNDArray< std::complex<float> >& ref = recon_obj.ref_calib_;

                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_parallel_imaging_lamda             : " << this->spirit_parallel_imaging_lamda.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_image_reg_lamda                    : " << this->spirit_image_reg_lamda.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_data_fidelity_lamda                : " << this->spirit_data_fidelity_lamda.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_nl_iter_max                        : " << this->spirit_nl_iter_max.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_nl_iter_thres                      : " << this->spirit_nl_iter_thres.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_name                           : " << this->spirit_reg_name.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_level                          : " << this->spirit_reg_level.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_keep_approx_coeff              : " << this->spirit_reg_keep_approx_coeff.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_keep_redundant_dimension_coeff : " << this->spirit_reg_keep_redundant_dimension_coeff.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_proximity_across_cha           : " << this->spirit_reg_proximity_across_cha.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_use_coil_sen_map               : " << this->spirit_reg_use_coil_sen_map.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_RO_weighting_ratio             : " << this->spirit_reg_RO_weighting_ratio.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_E1_weighting_ratio             : " << this->spirit_reg_E1_weighting_ratio.value());
                GDEBUG_CONDITION_STREAM(this->verbose.value(), "spirit_reg_N_weighting_ratio              : " << this->spirit_reg_N_weighting_ratio.value());

                size_t slc, s;

                for (slc = 0; slc < SLC; slc++)
                {
                    for (s = 0; s < S; s++)
                    {
                        std::stringstream os;
                        os << "encoding_" << e << "_s" << s << "_slc" << slc;
                        std::string suffix_2DT = os.str();

                        // ------------------------------

                        std::complex<float>* pKspace = &kspace(0, 0, 0, 0, 0, s, slc);
                        hoNDArray< std::complex<float> > kspace2DT(RO, E1, E2, dstCHA, N, 1, 1, pKspace);

                        // ------------------------------

                        long long kernelS = s;
                        if (kernelS >= (long long)ref_S) kernelS = (long long)ref_S - 1;

                        std::complex<float>* pKIm = &recon_obj.kernelIm2D_(0, 0, 0, 0, 0, kernelS, slc);
                        hoNDArray< std::complex<float> > kIm2DT(RO, E1, srcCHA, dstCHA, ref_N, 1, 1, pKIm);

                        // ------------------------------

                        std::complex<float>* pRef = &ref(0, 0, 0, 0, 0, kernelS, slc);
                        hoNDArray< std::complex<float> > ref2DT(ref.get_size(0), ref.get_size(1), ref.get_size(2), dstCHA, ref_N, 1, 1, pRef);

                        // ------------------------------

                        hoNDArray< std::complex<float> > coilMap2DT;
                        if (recon_obj.coil_map_.get_size(6) == SLC)
                        {
                            size_t coil_S = recon_obj.coil_map_.get_size(5);
                            std::complex<float>* pCoilMap = &recon_obj.coil_map_(0, 0, 0, 0, 0, ((s>=coil_S) ? coil_S-1 : s), slc);
                            coilMap2DT.create(RO, E1, E2, dstCHA, ref_N, 1, 1, pCoilMap);
                        }

                        // ------------------------------

                        std::complex<float>* pRes = &res(0, 0, 0, 0, 0, s, slc);
                        hoNDArray< std::complex<float> > res2DT(RO, E1, E2, dstCHA, N, 1, 1, pRes);

                        // ------------------------------

                        if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(kspace2DT, debug_folder_full_path_ + "kspace2DT_nl_spirit_" + suffix_2DT); }
                        if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(kIm2DT, debug_folder_full_path_ + "kIm2DT_nl_spirit_" + suffix_2DT); }
                        if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(ref2DT, debug_folder_full_path_ + "ref2DT_nl_spirit_" + suffix_2DT); }

                        // ------------------------------

                        std::string timing_str = "SPIRIT, Non-linear unwrapping, 2DT_" + suffix_2DT;
                        if (this->perform_timing.value()) timer.start(timing_str.c_str());
                        this->perform_nonlinear_spirit_unwrapping(kspace2DT, kIm2DT, ref2DT, coilMap2DT, res2DT, e);
                        if (this->perform_timing.value()) timer.stop();

                        if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(res2DT, debug_folder_full_path_ + "res_nl_spirit_2DT_" + suffix_2DT); }
                    }
                }
            }

            // ---------------------------------------------------------------------
            // compute coil combined images
            // ---------------------------------------------------------------------
            if (this->perform_timing.value()) timer.start("SPIRIT Non linear, coil combination ... ");
            this->perform_spirit_coil_combine(recon_obj);
            if (this->perform_timing.value()) timer.stop();

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(recon_obj.recon_res_.data_, debug_folder_full_path_ + "unwrappedIm_" + suffix); }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianNonLinearSpirit2DTGadget::perform_unwrapping(...) ... ");
        }
    }

    class solverCallBack : public hoGdSolverCallBack< hoNDArray< std::complex<float> >, hoWavelet2DTOperator< std::complex<float> > >
    {
        public:
        typedef hoGdSolverCallBack< hoNDArray< std::complex<float> >, hoWavelet2DTOperator< std::complex<float> > > BaseClass;

        solverCallBack() : BaseClass() {}
        virtual ~solverCallBack() {}

        void execute(const hoNDArray< std::complex<float> >& b, hoNDArray< std::complex<float> >& x)
        {
            typedef hoSPIRIT2DTDataFidelityOperator< std::complex<float> > SpiritOperType;
            SpiritOperType* pOper = dynamic_cast<SpiritOperType*> (this->solver_->oper_system_);
            pOper->restore_acquired_kspace(x);
        }
    };

    void GenericReconCartesianNonLinearSpirit2DTGadget::perform_nonlinear_spirit_unwrapping(hoNDArray< std::complex<float> >& kspace, 
        hoNDArray< std::complex<float> >& kerIm, hoNDArray< std::complex<float> >& ref2DT, hoNDArray< std::complex<float> >& coilMap2DT, hoNDArray< std::complex<float> >& res, size_t e)
    {
        try
        {
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

            hoNDArray< std::complex<float> > kspaceLinear(kspace);
            res = kspace;

            // detect whether random sampling is used
            bool use_random_sampling = false;
            std::vector<long long> sampled_step_size;
            long long n, e1;
            for (n=0; n<(long long)N; n++)
            {
                long long prev_sampled_line = -1;
                for (e1=0; e1<(long long)E1; e1++)
                {
                    if(std::abs(kspace(RO/2, e1, 0, 0, 0, 0, 0))>0 && std::abs(kspace(RO/2, e1, 0, CHA-1, 0, 0, 0))>0)
                    {
                        if(prev_sampled_line>0)
                        {
                            sampled_step_size.push_back(e1 - prev_sampled_line);
                        }

                        prev_sampled_line = e1;
                    }
                }
            }

            if(sampled_step_size.size()>4)
            {
                size_t s;
                for (s=2; s<sampled_step_size.size()-1; s++)
                {
                    if(sampled_step_size[s]!=sampled_step_size[s-1])
                    {
                        use_random_sampling = true;
                        break;
                    }
                }
            }

            if(use_random_sampling)
            {
                GDEBUG_STREAM("SPIRIT Non linear, random sampling is detected ... ");
            }

            Gadgetron::GadgetronTimer timer(false);

            boost::shared_ptr< hoNDArray< std::complex<float> > > coilMap;

            bool hasCoilMap = false;
            if (coilMap2DT.get_size(0) == RO && coilMap2DT.get_size(1) == E1 && coilMap2DT.get_size(3)==CHA)
            {
                if (ref_N < N)
                {
                    coilMap = boost::shared_ptr< hoNDArray< std::complex<float> > >(new hoNDArray< std::complex<float> >(RO, E1, CHA, coilMap2DT.begin()));
                }
                else
                {
                    coilMap = boost::shared_ptr< hoNDArray< std::complex<float> > >(new hoNDArray< std::complex<float> >(RO, E1, CHA, ref_N, coilMap2DT.begin()));
                }

                hasCoilMap = true;
            }

            hoNDArray<float> gFactor;
            float gfactorMedian = 0;

            float smallest_eigen_value(0);

            // -----------------------------------------------------
            // estimate gfactor
            // -----------------------------------------------------

            // mean over N
            hoNDArray< std::complex<float> > meanKSpace;

            if(calib_mode_[e]==ISMRMRD_interleaved)
            {
                Gadgetron::compute_averaged_data_N_S(kspace, true, true, true, meanKSpace);
            }
            else
            {
                Gadgetron::compute_averaged_data_N_S(ref2DT, true, true, true, meanKSpace);
            }

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(meanKSpace, debug_folder_full_path_ + "spirit_nl_2DT_meanKSpace"); }

            hoNDArray< std::complex<float> > acsSrc(meanKSpace.get_size(0), meanKSpace.get_size(1), CHA, meanKSpace.begin());
            hoNDArray< std::complex<float> > acsDst(meanKSpace.get_size(0), meanKSpace.get_size(1), CHA, meanKSpace.begin());

            double grappa_reg_lamda = 0.0005;
            size_t kRO = 5;
            size_t kE1 = 4;

            hoNDArray< std::complex<float> > convKer;
            hoNDArray< std::complex<float> > kIm(RO, E1, CHA, CHA);

            Gadgetron::grappa2d_calib_convolution_kernel(acsSrc, acsDst, (size_t)this->acceFactorE1_[e], grappa_reg_lamda, kRO, kE1, convKer);
            Gadgetron::grappa2d_image_domain_kernel(convKer, RO, E1, kIm);

            hoNDArray< std::complex<float> > unmixC;

            if(hasCoilMap)
            {
                Gadgetron::grappa2d_unmixing_coeff(kIm, *coilMap, (size_t)acceFactorE1_[e], unmixC, gFactor);

                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array(gFactor, debug_folder_full_path_ + "spirit_nl_2DT_gFactor");

                hoNDArray<float> gfactorSorted(gFactor);
                std::sort(gfactorSorted.begin(), gfactorSorted.begin()+RO*E1);
                gfactorMedian = gFactor((RO*E1 / 2));

                GDEBUG_STREAM("SPIRIT Non linear, the median gfactor is found to be : " << gfactorMedian);
            }

            if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(kIm, debug_folder_full_path_ + "spirit_nl_2DT_kIm");

            hoNDArray< std::complex<float> > complexIm;

            // compute linear solution as the initialization
            if(use_random_sampling)
            {
                if (this->perform_timing.value()) timer.start("SPIRIT Non linear, perform linear spirit recon ... ");
                this->perform_spirit_unwrapping(kspace, kerIm, kspaceLinear);
                if (this->perform_timing.value()) timer.stop();
            }
            else
            {
                if (this->perform_timing.value()) timer.start("SPIRIT Non linear, perform linear recon ... ");

                //size_t ref2DT_RO = ref2DT.get_size(0);
                //size_t ref2DT_E1 = ref2DT.get_size(1);

                //// mean over N
                //hoNDArray< std::complex<float> > meanKSpace;
                //Gadgetron::sum_over_dimension(ref2DT, meanKSpace, 4);

                //if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(meanKSpace, debug_folder_full_path_ + "spirit_nl_2DT_meanKSpace"); }

                //hoNDArray< std::complex<float> > acsSrc(ref2DT_RO, ref2DT_E1, CHA, meanKSpace.begin());
                //hoNDArray< std::complex<float> > acsDst(ref2DT_RO, ref2DT_E1, CHA, meanKSpace.begin());

                //double grappa_reg_lamda = 0.0005;
                //size_t kRO = 5;
                //size_t kE1 = 4;

                //hoNDArray< std::complex<float> > convKer;
                //hoNDArray< std::complex<float> > kIm(RO, E1, CHA, CHA);

                //Gadgetron::grappa2d_calib_convolution_kernel(acsSrc, acsDst, (size_t)this->acceFactorE1_[e], grappa_reg_lamda, kRO, kE1, convKer);
                //Gadgetron::grappa2d_image_domain_kernel(convKer, RO, E1, kIm);

                //if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(kIm, debug_folder_full_path_ + "spirit_nl_2DT_kIm");

                Gadgetron::hoNDFFT<float>::instance()->ifft2c(kspace, complex_im_recon_buf_);
                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(complex_im_recon_buf_, debug_folder_full_path_ + "spirit_nl_2DT_aliasedImage");

                hoNDArray< std::complex<float> > resKSpace(RO, E1, CHA, N);
                hoNDArray< std::complex<float> > aliasedImage(RO, E1, CHA, N, complex_im_recon_buf_.begin());
                Gadgetron::grappa2d_image_domain_unwrapping_aliased_image(aliasedImage, kIm, resKSpace);

                if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(resKSpace, debug_folder_full_path_ + "spirit_nl_2DT_linearImage");

                Gadgetron::hoNDFFT<float>::instance()->fft2c(resKSpace);

                memcpy(kspaceLinear.begin(), resKSpace.begin(), resKSpace.get_number_of_bytes());

                Gadgetron::apply_unmix_coeff_aliased_image(aliasedImage, unmixC, complexIm);

                if (this->perform_timing.value()) timer.stop();
            }

            if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(kspaceLinear, debug_folder_full_path_ + "spirit_nl_2DT_kspaceLinear");

            if(hasCoilMap)
            {
                if(N>=spirit_reg_minimal_num_images_for_noise_floor.value())
                {
                    // estimate the noise level

                    if(use_random_sampling)
                    {
                        Gadgetron::hoNDFFT<float>::instance()->ifft2c(kspaceLinear, complex_im_recon_buf_);

                        hoNDArray< std::complex<float> > complexLinearImage(RO, E1, CHA, N, complex_im_recon_buf_.begin());

                        Gadgetron::coil_combine(complexLinearImage, *coilMap, 2, complexIm);
                    }

                    if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(complexIm, debug_folder_full_path_ + "spirit_nl_2DT_linearImage_complexIm");

                    // if N is sufficiently large, we can estimate the noise floor by the smallest eigen value
                    hoMatrix< std::complex<float> > data;
                    data.createMatrix(RO*E1, N, complexIm.begin(), false);

                    hoNDArray< std::complex<float> > eigenVectors, eigenValues, eigenVectorsPruned;

                    // compute eigen
                    hoNDKLT< std::complex<float> > klt;
                    klt.prepare(data, (size_t)1, (size_t)0);
                    klt.eigen_value(eigenValues);

                    if (this->verbose.value())
                    {
                        GDEBUG_STREAM("SPIRIT Non linear, computes eigen values for all 2D kspaces ... ");
                        std::stringstream stream; 
                        eigenValues.print(stream);
                        GDEBUG(stream.str().c_str());

                        for (size_t i = 0; i<eigenValues.get_size(0); i++)
                        {
                            GDEBUG_STREAM(i << " = " << eigenValues(i));
                        }
                    }

                    smallest_eigen_value = std::sqrt( std::abs(eigenValues(N - 1).real()) / (RO*E1) );
                    GDEBUG_STREAM("SPIRIT Non linear, the smallest eigen value is : " << smallest_eigen_value);
                }
            }

            // perform nonlinear reconstruction
            {
                boost::shared_ptr<hoNDArray< std::complex<float> > > ker(new hoNDArray< std::complex<float> >(RO, E1, CHA, CHA, ref_N, kerIm.begin()));
                boost::shared_ptr<hoNDArray< std::complex<float> > > acq(new hoNDArray< std::complex<float> >(RO, E1, CHA, N, kspace.begin()));
                hoNDArray< std::complex<float> > kspaceInitial(RO, E1, CHA, N, kspaceLinear.begin());
                hoNDArray< std::complex<float> > res2DT(RO, E1, CHA, N, res.begin());

                if (this->spirit_data_fidelity_lamda.value() > 0)
                {
                    GDEBUG_STREAM("Start the NL SPIRIT data fidelity iteration - regularization strength : " << this->spirit_image_reg_lamda.value()
                                    << " - number of iteration : "                      << this->spirit_nl_iter_max.value()
                                    << " - proximity across cha : "                     << this->spirit_reg_proximity_across_cha.value()
                                    << " - redundant dimension weighting ratio : "      << this->spirit_reg_N_weighting_ratio.value()
                                    << " - using coil sen map : "                       << this->spirit_reg_use_coil_sen_map.value()
                                    << " - iter thres : "                               << this->spirit_nl_iter_thres.value()
                                    << " - wavelet name : "                             << this->spirit_reg_name.value()
                                    );

                    typedef hoGdSolver< hoNDArray< std::complex<float> >, hoWavelet2DTOperator< std::complex<float> > > SolverType;
                    SolverType solver;
                    solver.iterations_ = this->spirit_nl_iter_max.value();
                    solver.set_output_mode(this->spirit_print_iter.value() ? SolverType::OUTPUT_VERBOSE : SolverType::OUTPUT_SILENT);
                    solver.grad_thres_ = this->spirit_nl_iter_thres.value();

                    if(spirit_reg_estimate_noise_floor.value() && std::abs(smallest_eigen_value)>0)
                    {
                        solver.scale_factor_ = smallest_eigen_value;
                        solver.proximal_strength_ratio_ = this->spirit_image_reg_lamda.value() * gfactorMedian;

                        GDEBUG_STREAM("SPIRIT Non linear, eigen value is used to derive the regularization strength : " << solver.proximal_strength_ratio_ << " - smallest eigen value : " << solver.scale_factor_);
                    }
                    else
                    {
                        solver.proximal_strength_ratio_ = this->spirit_image_reg_lamda.value();
                    }

                    boost::shared_ptr< hoNDArray< std::complex<float> > > x0 = boost::make_shared< hoNDArray< std::complex<float> > >(kspaceInitial);
                    solver.set_x0(x0);

                    // parallel imaging term
                    std::vector<size_t> dims;
                    acq->get_dimensions(dims);
                    hoSPIRIT2DTDataFidelityOperator< std::complex<float> > spirit(dims);
                    spirit.set_forward_kernel(*ker, false);
                    spirit.set_acquired_points(*acq);

                    // image reg term
                    hoWavelet2DTOperator< std::complex<float> > wav3DOperator(dims);
                    wav3DOperator.set_acquired_points(*acq);
                    wav3DOperator.scale_factor_first_dimension_ = this->spirit_reg_RO_weighting_ratio.value();
                    wav3DOperator.scale_factor_second_dimension_ = this->spirit_reg_E1_weighting_ratio.value();
                    wav3DOperator.scale_factor_third_dimension_ = this->spirit_reg_N_weighting_ratio.value();
                    wav3DOperator.with_approx_coeff_ = !this->spirit_reg_keep_approx_coeff.value();
                    wav3DOperator.change_coeffcients_third_dimension_boundary_ = !this->spirit_reg_keep_redundant_dimension_coeff.value();
                    wav3DOperator.proximity_across_cha_ = this->spirit_reg_proximity_across_cha.value();
                    wav3DOperator.no_null_space_ = true;
                    wav3DOperator.input_in_kspace_ = true;
                    wav3DOperator.select_wavelet(this->spirit_reg_name.value());

                    if (this->spirit_reg_use_coil_sen_map.value() && hasCoilMap)
                    {
                        wav3DOperator.coil_map_ = *coilMap;
                    }

                    // set operators

                    solver.oper_system_ = &spirit;
                    solver.oper_reg_ = &wav3DOperator;

                    if (this->perform_timing.value()) timer.start("NonLinear SPIRIT solver for 2DT with data fidelity ... ");
                    solver.solve(*acq, res2DT);
                    if (this->perform_timing.value()) timer.stop();

                    if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(res2DT, debug_folder_full_path_ + "spirit_nl_2DT_data_fidelity_res");
                }
                else
                {
                    GDEBUG_STREAM("Start the NL SPIRIT iteration with regularization strength : "<< this->spirit_image_reg_lamda.value()
                                    << " - number of iteration : " << this->spirit_nl_iter_max.value()
                                    << " - proximity across cha : " << this->spirit_reg_proximity_across_cha.value()
                                    << " - redundant dimension weighting ratio : " << this->spirit_reg_N_weighting_ratio.value()
                                    << " - using coil sen map : " << this->spirit_reg_use_coil_sen_map.value()
                                    << " - iter thres : " << this->spirit_nl_iter_thres.value()
                                    << " - wavelet name : " << this->spirit_reg_name.value()
                                    );

                    typedef hoGdSolver< hoNDArray< std::complex<float> >, hoWavelet2DTOperator< std::complex<float> > > SolverType;
                    SolverType solver;
                    solver.iterations_ = this->spirit_nl_iter_max.value();
                    solver.set_output_mode(this->spirit_print_iter.value() ? SolverType::OUTPUT_VERBOSE : SolverType::OUTPUT_SILENT);
                    solver.grad_thres_ = this->spirit_nl_iter_thres.value();

                    if(spirit_reg_estimate_noise_floor.value() && std::abs(smallest_eigen_value)>0)
                    {
                        solver.scale_factor_ = smallest_eigen_value;
                        solver.proximal_strength_ratio_ = this->spirit_image_reg_lamda.value() * gfactorMedian;

                        GDEBUG_STREAM("SPIRIT Non linear, eigen value is used to derive the regularization strength : " << solver.proximal_strength_ratio_ << " - smallest eigen value : " << solver.scale_factor_);
                    }
                    else
                    {
                        solver.proximal_strength_ratio_ = this->spirit_image_reg_lamda.value();
                    }

                    boost::shared_ptr< hoNDArray< std::complex<float> > > x0 = boost::make_shared< hoNDArray< std::complex<float> > >(kspaceInitial);
                    solver.set_x0(x0);

                    // parallel imaging term
                    std::vector<size_t> dims;
                    acq->get_dimensions(dims);

                    hoSPIRIT2DTOperator< std::complex<float> > spirit(dims);
                    spirit.set_forward_kernel(*ker, false);
                    spirit.set_acquired_points(*acq);
                    spirit.no_null_space_ = true;
                    spirit.use_non_centered_fft_ = false;

                    // image reg term
                    std::vector<size_t> dim;
                    acq->get_dimensions(dim);

                    hoWavelet2DTOperator< std::complex<float> > wav3DOperator(dim);
                    wav3DOperator.set_acquired_points(*acq);
                    wav3DOperator.scale_factor_first_dimension_ = this->spirit_reg_RO_weighting_ratio.value();
                    wav3DOperator.scale_factor_second_dimension_ = this->spirit_reg_E1_weighting_ratio.value();
                    wav3DOperator.scale_factor_third_dimension_ = this->spirit_reg_N_weighting_ratio.value();
                    wav3DOperator.with_approx_coeff_ = !this->spirit_reg_keep_approx_coeff.value();
                    wav3DOperator.change_coeffcients_third_dimension_boundary_ = !this->spirit_reg_keep_redundant_dimension_coeff.value();
                    wav3DOperator.proximity_across_cha_ = this->spirit_reg_proximity_across_cha.value();
                    wav3DOperator.no_null_space_ = true;
                    wav3DOperator.input_in_kspace_ = true;
                    wav3DOperator.select_wavelet(this->spirit_reg_name.value());

                    if (this->spirit_reg_use_coil_sen_map.value() && hasCoilMap)
                    {
                        wav3DOperator.coil_map_ = *coilMap;
                    }

                    // set operators
                    solver.oper_system_ = &spirit;
                    solver.oper_reg_ = &wav3DOperator;

                    // set call back
                    solverCallBack cb;
                    cb.solver_ = &solver;
                    solver.call_back_ = &cb;

                    hoNDArray< std::complex<float> > b(kspaceInitial);
                    Gadgetron::clear(b);

                    if (this->perform_timing.value()) timer.start("NonLinear SPIRIT solver for 2DT ... ");
                    solver.solve(b, res2DT);
                    if (this->perform_timing.value()) timer.stop();

                    if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(res2DT, debug_folder_full_path_ + "spirit_nl_2DT_res");

                    spirit.restore_acquired_kspace(kspace, res2DT);

                    if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(res2DT, debug_folder_full_path_ + "spirit_nl_2DT_res_restored");
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianNonLinearSpirit2DTGadget::perform_nonlinear_spirit_unwrapping(...) ... ");
        }
    }

    GADGET_FACTORY_DECLARE(GenericReconCartesianNonLinearSpirit2DTGadget)
}
