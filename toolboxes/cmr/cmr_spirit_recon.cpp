/** \file   cmr_spirit_recon.cpp
    \brief  Implement some functionalities commonly used in cmr applications for spirit recon
    \author Hui Xue
*/

#include "cmr_spirit_recon.h"

#include "mri_core_utility.h"
#include "mri_core_spirit.h"
#include "hoNDArray_reductions.h"
#include "hoNDFFT.h"
#include "hoSPIRIT2DOperator.h"
#include "hoLsqrSolver.h"
#include "hoSPIRIT2DTDataFidelityOperator.h"
#include "hoWavelet2DTOperator.h"
#include "hoGdSolver.h"
#include <boost/make_shared.hpp>

namespace Gadgetron {

    template <typename T>
    void perform_spirit_recon_linear_2DT(const Gadgetron::hoNDArray<T>& kspace, size_t startE1, size_t endE1, const Gadgetron::hoNDArray<T>& kerIm,
                                const Gadgetron::hoNDArray<T>& kspaceInitial, Gadgetron::hoNDArray<T>& res, size_t iter_max, double iter_thres, bool print_iter)
    {
        try
        {
            size_t RO = kspace.get_size(0);
            size_t E1 = kspace.get_size(1);
            size_t CHA = kspace.get_size(2);
            size_t N = kspace.get_size(3);
            size_t S = kspace.get_size(4);

            size_t ref_N = kerIm.get_size(4);
            size_t ref_S = kerIm.get_size(5);

            if(startE1>=E1) startE1 = 0;
            if(endE1>=E1) startE1 = 0;
            if(startE1>endE1)
            {
                startE1 = 0;
                endE1 = E1 - 1;
            }

            long long num = N*S;

            hoNDArray<T> ker_Shifted(kerIm);
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(kerIm, ker_Shifted);

            hoNDArray<T> kspace_Shifted;
            kspace_Shifted = kspace;
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(kspace, kspace_Shifted);

            hoNDArray<T> kspace_initial_Shifted;
            bool hasInitial = false;
            if ( kspaceInitial.dimensions_equal(kspace) )
            {
                kspace_initial_Shifted = kspaceInitial;
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(kspaceInitial, kspace_initial_Shifted);
                hasInitial = true;
            }

            res = kspace;

#ifdef USE_OMP
            int numThreads = (int)num;
            if (numThreads > omp_get_num_procs()) numThreads = omp_get_num_procs();
            GDEBUG_CONDITION_STREAM(print_iter, "numThreads : " << numThreads);
#endif // USE_OMP

            long long ii;

            std::vector<size_t> dim(3, 1);
            dim[0] = RO;
            dim[1] = E1;
            dim[2] = CHA;

#pragma omp parallel default(none) private(ii) shared(num, N, S, RO, E1, CHA, dim, startE1, endE1, ref_N, ref_S, kspace, res, kspace_Shifted, ker_Shifted, kspace_initial_Shifted, hasInitial, iter_max, iter_thres, print_iter) num_threads(numThreads) if(num>1)
            {
                boost::shared_ptr< hoSPIRIT2DOperator< T > > oper(new hoSPIRIT2DOperator< T >(dim));
                hoSPIRIT2DOperator< T >& spirit = *oper;
                spirit.use_non_centered_fft_ = true;
                spirit.no_null_space_ = false;

                if (ref_N == 1 && ref_S == 1)
                {
                    boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray< T >(RO, E1, CHA, CHA, ker_Shifted.begin()));
                    spirit.set_forward_kernel(*ker, false);
                }

                hoLsqrSolver< T > cgSolver;
                cgSolver.set_tc_tolerance((float)iter_thres);
                cgSolver.set_max_iterations( (unsigned int)iter_max);
                cgSolver.set_output_mode(print_iter ? hoLsqrSolver< T >::OUTPUT_VERBOSE : hoLsqrSolver< T >::OUTPUT_SILENT);
                cgSolver.set_encoding_operator(oper);

                hoNDArray< T > b(RO, E1, CHA);
                hoNDArray< T > unwarppedKSpace(RO, E1, CHA);

#pragma omp for
                for (ii = 0; ii < num; ii++)
                {
                    size_t s = ii / N;
                    size_t n = ii - s*N;

                    // check whether the kspace is undersampled
                    bool undersampled = false;
                    for (size_t e1 = startE1; e1 <= endE1; e1++)
                    {
                        if ((std::abs(kspace(RO / 2, e1, CHA - 1, n, s)) == 0)
                            && (std::abs(kspace(RO / 2, e1, 0, n, s)) == 0))
                        {
                            undersampled = true;
                            break;
                        }
                    }

                    T* pKpaceShifted = &(kspace_Shifted(0, 0, 0, n, s));
                    T* pRes = &(res(0, 0, 0, n, s));

                    if (!undersampled)
                    {
                        memcpy(pRes, pKpaceShifted, sizeof(T)*RO*E1*CHA);
                        continue;
                    }

                    long long kernelN = n;
                    if (kernelN >= (long long)ref_N) kernelN = (long long)ref_N - 1;

                    long long kernelS = s;
                    if (kernelS >= (long long)ref_S) kernelS = (long long)ref_S - 1;

                    boost::shared_ptr< hoNDArray< T > > acq(new hoNDArray< T >(RO, E1, CHA, pKpaceShifted));
                    spirit.set_acquired_points(*acq);

                    boost::shared_ptr< hoNDArray<T> > initialAcq;
                    if ( hasInitial )
                    {
                        initialAcq = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>(RO, E1, CHA, &kspace_initial_Shifted(0, 0, 0, n, s)));
                        cgSolver.set_x0(initialAcq);
                    }
                    else
                    {
                        cgSolver.set_x0(acq);
                    }

                    if (ref_N == 1 && ref_S == 1)
                    {
                        spirit.compute_righ_hand_side(*acq, b);
                        cgSolver.solve(&unwarppedKSpace, &b);
                    }
                    else
                    {
                        T* pKer = &(ker_Shifted(0, 0, 0, kernelN, kernelS));
                        boost::shared_ptr<hoNDArray< T > > ker(new hoNDArray< T >(RO, E1, CHA, CHA, pKer));
                        spirit.set_forward_kernel(*ker, false);

                        spirit.compute_righ_hand_side(*acq, b);
                        cgSolver.solve(&unwarppedKSpace, &b);
                    }

                    // restore the acquired points
                    spirit.restore_acquired_kspace(*acq, unwarppedKSpace);
                    memcpy(pRes, unwarppedKSpace.begin(), unwarppedKSpace.get_number_of_bytes());
                }
            }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fftshift2D(res, kspace_Shifted);
            res = kspace_Shifted;
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in perform_spirit_recon_linear_2DT(...) ... ");
        }
    }

    template void perform_spirit_recon_linear_2DT(const Gadgetron::hoNDArray< std::complex<float> >& kspace, size_t startE1, size_t endE1, const Gadgetron::hoNDArray< std::complex<float> >& kerIm, const Gadgetron::hoNDArray< std::complex<float> >& kspaceInitial, Gadgetron::hoNDArray< std::complex<float> >& res, size_t iter_max, double iter_thres, bool print_iter);
    template void perform_spirit_recon_linear_2DT(const Gadgetron::hoNDArray< std::complex<double> >& kspace, size_t startE1, size_t endE1, const Gadgetron::hoNDArray< std::complex<double> >& kerIm, const Gadgetron::hoNDArray< std::complex<double> >& kspaceInitial, Gadgetron::hoNDArray< std::complex<double> >& res, size_t iter_max, double iter_thres, bool print_iter);

    // ---------------------------------------------------------------------

    template <typename T>
    class sCB : public hoGdSolverCallBack< hoNDArray< T >, hoWavelet2DTOperator< T > >
    {
    public:
        typedef hoGdSolverCallBack< hoNDArray<T>, hoWavelet2DTOperator<T> > BaseClass;

        sCB() : BaseClass() {}
        virtual ~sCB() {}

        void execute(const hoNDArray<T>& b, hoNDArray<T>& x)
        {
            typedef hoSPIRIT2DTDataFidelityOperator<T> SpiritOperType;
            SpiritOperType* pOper = dynamic_cast<SpiritOperType*> (this->solver_->oper_system_);
            pOper->restore_acquired_kspace(x);
        }
    };

    template <typename T>
    void perform_spirit_recon_non_linear_2DT(const Gadgetron::hoNDArray<T>& kspace, const Gadgetron::hoNDArray<T>& kerIm,
                                            const Gadgetron::hoNDArray<T>& coil_map, const Gadgetron::hoNDArray<T>& kspaceLinear, Gadgetron::hoNDArray<T>& res,
                                            size_t iter_max, double iter_thres, double data_fidelity_lamda, double image_reg_lamda, double reg_N_weighting_ratio,
                                            bool reg_use_coil_sen_map, bool reg_with_approx_coeff, const std::string& wav_name, bool print_iter)
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t CHA = kspace.get_size(2);
        size_t N = kspace.get_size(3);
        size_t S = kspace.get_size(4);

        size_t ref_N = kerIm.get_size(4);
        size_t ref_S = kerIm.get_size(5);

        res = kspace;

        size_t s;

        for (s=0; s<S; s++)
        {
            boost::shared_ptr< hoNDArray< T > > coilMap;

            bool hasCoilMap = false;
            if (coil_map.get_size(0) == RO && coil_map.get_size(1) == E1 && coil_map.get_size(2)==CHA)
            {
                if (ref_N < N)
                {
                    coilMap = boost::shared_ptr< hoNDArray< T > >(new hoNDArray< T >(RO, E1, CHA, const_cast<T*>(coil_map.begin()) ));
                }
                else
                {
                    coilMap = boost::shared_ptr< hoNDArray< T > >(new hoNDArray< T >(RO, E1, CHA, ref_N, const_cast<T*>(coil_map.begin()) ));
                }

                hasCoilMap = true;
            }

            size_t s_used = s;
            if(s_used>=ref_S) s_used = ref_S;

            boost::shared_ptr<hoNDArray< T > > ker(new hoNDArray< T >(RO, E1, CHA, CHA, ref_N, const_cast<T*>(kerIm.begin())+s_used*RO*E1*CHA*CHA*ref_N) );
            boost::shared_ptr<hoNDArray< T > > acq(new hoNDArray< T >(RO, E1, CHA, N, const_cast<T*>(kspace.begin())+s*RO*E1*CHA*N) );
            hoNDArray< T > kspaceInitial(RO, E1, CHA, N, const_cast<T*>(kspaceLinear.begin())+s*RO*E1*CHA*N );
            hoNDArray< T > res2DT(RO, E1, CHA, N, res.begin()+s*RO*E1*CHA*N );

            if (data_fidelity_lamda > 0)
            {
                GDEBUG_STREAM("Start the NonLinear SPIRIT data fidelity iteration - regularization strength : "
                    << image_reg_lamda
                    << " - number of iteration : "                      << iter_max
                    << " - proximity across cha : "                     << false
                    << " - redundant dimension weighting ratio : "      << reg_N_weighting_ratio
                    << " - using coil sen map : "                       << reg_use_coil_sen_map
                    << " - with approx coeff : "                        << reg_with_approx_coeff
                    << " - iter thres : "                               << iter_thres
                    << " - wavelet name : "                             << wav_name
                    );

                typedef hoGdSolver< hoNDArray< T >, hoWavelet2DTOperator< T > > SolverType;
                SolverType solver;
                solver.iterations_ = iter_max;
                solver.set_output_mode(print_iter ? SolverType::OUTPUT_VERBOSE : SolverType::OUTPUT_SILENT);
                solver.grad_thres_ = iter_thres;
                solver.proximal_strength_ratio_ = image_reg_lamda;

                boost::shared_ptr< hoNDArray< T > > x0 = boost::make_shared< hoNDArray< T > >(kspaceInitial);
                solver.set_x0(x0);

                // parallel imaging term
                std::vector<size_t> dims;
                acq->get_dimensions(dims);
                hoSPIRIT2DTDataFidelityOperator< T > spirit(dims);
                spirit.set_forward_kernel(*ker, false);
                spirit.set_acquired_points(*acq);

                // image reg term
                hoWavelet2DTOperator< T > wav3DOperator(dims);
                wav3DOperator.set_acquired_points(*acq);
                wav3DOperator.scale_factor_first_dimension_ = 1;
                wav3DOperator.scale_factor_second_dimension_ = 1;
                wav3DOperator.scale_factor_third_dimension_ = reg_N_weighting_ratio;
                wav3DOperator.with_approx_coeff_ = reg_with_approx_coeff;
                wav3DOperator.change_coeffcients_third_dimension_boundary_ = true;
                wav3DOperator.proximity_across_cha_ = false;
                wav3DOperator.no_null_space_ = true;
                wav3DOperator.input_in_kspace_ = true;
                wav3DOperator.select_wavelet(wav_name);

                if (reg_use_coil_sen_map && hasCoilMap)
                {
                    wav3DOperator.coil_map_ = *coilMap;
                }

                // set operators

                solver.oper_system_ = &spirit;
                solver.oper_reg_ = &wav3DOperator;

                solver.solve(*acq, res2DT);
            }
            else
            {
                GDEBUG_STREAM("Start the NonLinear SPIRIT iteration with regularization strength : "
                    << image_reg_lamda
                    << " - number of iteration : " << iter_max
                    << " - proximity across cha : " << false
                    << " - redundant dimension weighting ratio : " << reg_N_weighting_ratio
                    << " - using coil sen map : " << reg_use_coil_sen_map
                    << " - with approx coeff : "  << reg_with_approx_coeff
                    << " - iter thres : " << iter_thres
                    << " - wavelet name : " << wav_name
                    );

                typedef hoGdSolver< hoNDArray< T >, hoWavelet2DTOperator< T > > SolverType;
                SolverType solver;
                solver.iterations_ = iter_max;
                solver.set_output_mode(print_iter ? SolverType::OUTPUT_VERBOSE : SolverType::OUTPUT_SILENT);
                solver.grad_thres_ = iter_thres;
                solver.proximal_strength_ratio_ = image_reg_lamda;

                boost::shared_ptr< hoNDArray< T > > x0 = boost::make_shared< hoNDArray< T > >(kspaceInitial);
                solver.set_x0(x0);

                // parallel imaging term
                std::vector<size_t> dims;
                acq->get_dimensions(dims);

                hoSPIRIT2DTOperator< T > spirit(dims);
                spirit.set_forward_kernel(*ker, false);
                spirit.set_acquired_points(*acq);
                spirit.no_null_space_ = true;
                spirit.use_non_centered_fft_ = false;

                // image reg term
                std::vector<size_t> dim;
                acq->get_dimensions(dim);

                hoWavelet2DTOperator< T > wav3DOperator(dim);
                wav3DOperator.set_acquired_points(*acq);
                wav3DOperator.scale_factor_first_dimension_ = 1;
                wav3DOperator.scale_factor_second_dimension_ = 1;
                wav3DOperator.scale_factor_third_dimension_ = reg_N_weighting_ratio;
                wav3DOperator.with_approx_coeff_ = reg_with_approx_coeff;
                wav3DOperator.change_coeffcients_third_dimension_boundary_ = true;
                wav3DOperator.proximity_across_cha_ = false;
                wav3DOperator.no_null_space_ = true;
                wav3DOperator.input_in_kspace_ = true;
                wav3DOperator.select_wavelet(wav_name);

                if (reg_use_coil_sen_map && hasCoilMap)
                {
                    wav3DOperator.coil_map_ = *coilMap;
                }

                // set operators
                solver.oper_system_ = &spirit;
                solver.oper_reg_ = &wav3DOperator;

                // set call back
                sCB<T> cb;
                cb.solver_ = &solver;
                solver.call_back_ = &cb;

                hoNDArray< T > b(kspaceInitial);
                Gadgetron::clear(b);

                solver.solve(b, res2DT);
                // if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(res2DT, debug_folder_full_path_ + "spirit_nl_2DT_res");

                spirit.restore_acquired_kspace(*acq, res2DT);

                // if (!debug_folder_full_path_.empty()) gt_exporter_.export_array_complex(res2DT, debug_folder_full_path_ + "spirit_nl_2DT_res_restored");
            }
        }
    }

    template void perform_spirit_recon_non_linear_2DT(const Gadgetron::hoNDArray< std::complex<float> >& kspace, const Gadgetron::hoNDArray< std::complex<float> >& kerIm,
        const Gadgetron::hoNDArray< std::complex<float> >& coilMap, const Gadgetron::hoNDArray< std::complex<float> >& kspaceInitial, Gadgetron::hoNDArray< std::complex<float> >& res,
        size_t iter_max, double iter_thres, double data_fidelity_lamda, double image_reg_lamda, double reg_N_weighting_ratio,
        bool reg_use_coil_sen_map, bool reg_with_approx_coeff, const std::string& wav_name, bool print_iter);

    template void perform_spirit_recon_non_linear_2DT(const Gadgetron::hoNDArray< std::complex<double> >& kspace, const Gadgetron::hoNDArray< std::complex<double> >& kerIm,
        const Gadgetron::hoNDArray< std::complex<double> >& coilMap, const Gadgetron::hoNDArray< std::complex<double> >& kspaceInitial, Gadgetron::hoNDArray< std::complex<double> >& res,
        size_t iter_max, double iter_thres, double data_fidelity_lamda, double image_reg_lamda, double reg_N_weighting_ratio,
        bool reg_use_coil_sen_map, bool reg_with_approx_coeff, const std::string& wav_name, bool print_iter);
}
