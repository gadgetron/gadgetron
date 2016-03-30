/** \file   gtPlusISMRMRDReconWorker2DTL1SPIRITNCG.h
    \brief  Implement the 2DT non-linear SPIRIT reconstruction using the non-linear CG solver
    \author Hui Xue
*/

#pragma once

#include "gtPlusISMRMRDReconWorker2DTSPIRIT.h"
#include "gtPlusNCGSolver.h"
#include "gtPlusDataFidelityOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorker2DTL1SPIRITNCG : public gtPlusReconWorker2DTSPIRIT<T>
{
public:

    typedef gtPlusReconWorker2DTSPIRIT<T> BaseClass;
    typedef typename realType<T>::Type value_type;

    gtPlusReconWorker2DTL1SPIRITNCG() : BaseClass() {}
    virtual ~gtPlusReconWorker2DTL1SPIRITNCG() {}

    virtual bool performUnwarppingImpl(gtPlusReconWorkOrder<T>* workOrder2DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res, size_t s);
    virtual bool performUnwarppingImpl(gtPlusReconJob2DT<T>& job);
    // virtual bool performUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& data);

    virtual bool autoReconParameter(gtPlusReconWorkOrder<T>* workOrder);

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::verbose_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;

    using BaseClass::buffer2DT_;
    using BaseClass::buffer2DT_unwrapping_;
    using BaseClass::buffer2DT_partial_fourier_;
    using BaseClass::buffer2DT_partial_fourier_kspaceIter_;
    using BaseClass::ref_src_;
    using BaseClass::ref_dst_;
    using BaseClass::data_dst_;
    using BaseClass::ref_coil_map_dst_;
    using BaseClass::startE1_;
    using BaseClass::endE1_;

    using BaseClass::spirit_;
};

template <typename T> 
bool gtPlusReconWorker2DTL1SPIRITNCG<T>::autoReconParameter(gtPlusReconWorkOrder<T>* workOrder)
{
    BaseClass::autoReconParameter(workOrder);

    gtPlusReconWorkOrder2DT<T>* workOrder2DT = dynamic_cast<gtPlusReconWorkOrder2DT<T>*>(workOrder);
    if ( workOrder2DT == NULL ) return false;

    if ( workOrder2DT->spirit_perform_linear_ )
    {
        if ( workOrder2DT->spirit_solve_symmetric_ )
        {
            workOrder2DT->spirit_image_reg_lamda_ = 0.0025;
            workOrder2DT->spirit_ncg_iter_thres_ = 0.0001;
        }
        else
        {
            workOrder2DT->spirit_image_reg_lamda_ = 0.0025;
            workOrder2DT->spirit_ncg_iter_thres_ = 0.0001;
        }
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DTL1SPIRITNCG<T>::
performUnwarppingImpl(gtPlusReconWorkOrder<T>* workOrder2DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res, size_t s)
{
    try
    {
        hoNDArray<T> kspaceLinear(kspace);
        res = kspace;

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace"); }

        bool performLinear = workOrder2DT->spirit_perform_linear_;
        if ( !workOrder2DT->spirit_perform_nonlinear_ ) performLinear = true;

        if ( performLinear )
        {
            if ( performTiming_ ) { gt_timer3_.start("NCG spirit linear solver for 2DT ... "); }
            GADGET_CHECK_RETURN_FALSE(BaseClass::performUnwarppingImpl(workOrder2DT, kspace, adj_forward_G_I, kspaceLinear, s));
            if ( performTiming_ ) { gt_timer3_.stop(); }
        }
        else
        {
            if ( workOrder2DT->kspace_initial_.get_number_of_elements() == kspace.get_number_of_elements() )
            {
                GDEBUG_STREAM("Start the iteration with the input initial kspace ... ");
                memcpy(kspaceLinear.begin(), workOrder2DT->kspace_initial_.begin(), kspace.get_number_of_bytes());
            }
            else
            {
                GDEBUG_STREAM("Start the iteration with the input kspace ... ");
            }
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceLinear, debugFolder_+"kspaceLinear"); }

        if ( workOrder2DT->spirit_perform_nonlinear_ )
        {
            size_t refN = adj_forward_G_I.get_size(4);

            size_t RO = kspace.get_size(0);
            size_t E1 = kspace.get_size(1);
            size_t N = kspace.get_size(3);

            size_t srcCHA = adj_forward_G_I.get_size(2);
            size_t dstCHA = adj_forward_G_I.get_size(3);

            if ( workOrder2DT->spirit_2D_scale_per_chunk_ )
            {
                typename realType<T>::Type scaleFactor = 1.0;
                Gadgetron::norm2(kspace, scaleFactor);
                scaleFactor /= (typename realType<T>::Type)( (RO*std::sqrt(double(srcCHA))) );

                workOrder2DT->spirit_ncg_scale_factor_ = scaleFactor;
            }

            // apply the scale
            if ( workOrder2DT->spirit_ncg_scale_factor_ > 0 )
            {
                Gadgetron::scal( static_cast<value_type>(1.0/workOrder2DT->spirit_ncg_scale_factor_), kspaceLinear);
                Gadgetron::scal( static_cast<value_type>(1.0/workOrder2DT->spirit_ncg_scale_factor_), kspace);
            }

            boost::shared_ptr< hoNDArray<T> > coilMapS;
            
            if ( workOrder2DT->coilMap_ )
            {
                if ( refN < N )
                {
                    coilMapS = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>(RO, E1, dstCHA, workOrder2DT->coilMap_->begin()));
                }
                else
                {
                    coilMapS = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>(RO, E1, dstCHA, refN, workOrder2DT->coilMap_->begin()+s*RO*E1*dstCHA*refN));
                }
            }

            if ( N > 1 )
            {
                // 2D+T
                boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(RO, E1, srcCHA, dstCHA, refN, adj_forward_G_I.begin()));
                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(RO, E1, srcCHA, N, kspace.begin()));

                gtPlusNCGSolver<hoNDArray<T>, hoNDArray<T>, gtPlusOperator<T> > ncgsolver;
                ncgsolver.iterMax_ = workOrder2DT->spirit_ncg_iter_max_;
                ncgsolver.printIter_ = workOrder2DT->spirit_ncg_print_iter_;
                ncgsolver.secantRatio_ = 1;
                ncgsolver.x0_ = &kspaceLinear;

                hoNDArray<T> b;

                if ( workOrder2DT->spirit_data_fidelity_lamda_ <= 0 )
                {
                    // parallel imaging term
                    gtPlusSPIRIT2DTOperator<T> spirit;
                    spirit.use_symmetric_spirit_ = false;
                    spirit.setForwardKernel(ker, true);
                    spirit.setAcquiredPoints(acq);

                    // L1 term
                    gtPlusWavelet3DOperator<T> wavNullSpace3DOperator;
                    wavNullSpace3DOperator.setAcquiredPoints(acq);
                    wavNullSpace3DOperator.scale_factor_first_dimension_ = (value_type)workOrder2DT->spirit_RO_enhancement_ratio_;
                    wavNullSpace3DOperator.scale_factor_second_dimension_ = (value_type)workOrder2DT->spirit_E1_enhancement_ratio_;
                    wavNullSpace3DOperator.scale_factor_third_dimension_ = (value_type)workOrder2DT->spirit_temporal_enhancement_ratio_;

                    if ( workOrder2DT->spirit_use_coil_sen_map_ && workOrder2DT->coilMap_ )
                    {
                        wavNullSpace3DOperator.setCoilSenMap(coilMapS);
                    }

                    // set operators
                    ncgsolver.add(spirit, T( (value_type)workOrder2DT->spirit_parallel_imaging_lamda_ ) );
                    ncgsolver.add(wavNullSpace3DOperator, T( (value_type)workOrder2DT->spirit_image_reg_lamda_ ) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 2DT ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_2DT_res"); }

                    spirit.restoreAcquiredKSpace(kspace, res);

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_2DT_res_restored"); }
                }
                else
                {
                    gtPlusSPIRITNoNullSpace2DTOperator<T> spirit_noNullSpace;
                    spirit_noNullSpace.use_symmetric_spirit_ = false;
                    spirit_noNullSpace.setForwardKernel(ker, true);
                    spirit_noNullSpace.setAcquiredPoints(acq);

                    gtPlusDataFidelityOperator<T> dataOper;
                    dataOper.setAcquiredPoints(acq);

                    gtPlusWaveletNoNullSpace3DOperator<T> wavNoNullSpace3DOperator;
                    wavNoNullSpace3DOperator.setAcquiredPoints(acq);
                    wavNoNullSpace3DOperator.scale_factor_first_dimension_ = (value_type)workOrder2DT->spirit_RO_enhancement_ratio_;
                    wavNoNullSpace3DOperator.scale_factor_second_dimension_ = (value_type)workOrder2DT->spirit_E1_enhancement_ratio_;
                    wavNoNullSpace3DOperator.scale_factor_third_dimension_ = (value_type)workOrder2DT->spirit_temporal_enhancement_ratio_;

                    if ( workOrder2DT->spirit_use_coil_sen_map_ && workOrder2DT->coilMap_ )
                    {
                        wavNoNullSpace3DOperator.setCoilSenMap(coilMapS);
                    }

                    ncgsolver.add(spirit_noNullSpace, T( (value_type)workOrder2DT->spirit_parallel_imaging_lamda_ ) );
                    ncgsolver.add(wavNoNullSpace3DOperator, T( (value_type)workOrder2DT->spirit_image_reg_lamda_ ) );
                    ncgsolver.add(dataOper, T( (value_type)workOrder2DT->spirit_data_fidelity_lamda_ ) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 2DT without null space ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_2DT_res_noNullSpace"); }
                }
            }
            else
            {
                // 2D
                boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(RO, E1, srcCHA, dstCHA, adj_forward_G_I.begin()));
                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(RO, E1, srcCHA, kspace.begin()));

                gtPlusNCGSolver<hoNDArray<T>, hoNDArray<T>, gtPlusOperator<T> > ncgsolver;
                ncgsolver.iterMax_ = workOrder2DT->spirit_ncg_iter_max_;
                ncgsolver.printIter_ = workOrder2DT->spirit_ncg_print_iter_;
                ncgsolver.secantRatio_ = 1;
                ncgsolver.x0_ = &kspaceLinear;

                hoNDArray<T> b;

                if ( workOrder2DT->spirit_data_fidelity_lamda_ <= 0 )
                {
                    // parallel imaging term
                    gtPlusSPIRIT2DOperator<T> spirit;
                    spirit.use_symmetric_spirit_ = false;
                    spirit.setForwardKernel(ker, true);
                    spirit.setAcquiredPoints(acq);

                    // L1 term
                    gtPlusWavelet2DOperator<T> wavNullSpace2DOperator;
                    wavNullSpace2DOperator.setAcquiredPoints(acq);

                    if ( workOrder2DT->spirit_use_coil_sen_map_ && workOrder2DT->coilMap_ )
                    {
                        wavNullSpace2DOperator.setCoilSenMap(coilMapS);
                    }

                    // set operators
                    ncgsolver.add(spirit, T( (value_type)workOrder2DT->spirit_parallel_imaging_lamda_ ) );
                    ncgsolver.add(wavNullSpace2DOperator, T( (value_type)workOrder2DT->spirit_image_reg_lamda_ ) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 2D ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_2D_res"); }

                    spirit.restoreAcquiredKSpace(kspace, res);

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, "ncg_spirit_2D_res_restored"); }
                }
                else
                {
                    gtPlusSPIRITNoNullSpace2DOperator<T> spirit_noNullSpace;
                    spirit_noNullSpace.use_symmetric_spirit_ = false;
                    spirit_noNullSpace.setForwardKernel(ker, true);
                    spirit_noNullSpace.setAcquiredPoints(acq);

                    gtPlusDataFidelityOperator<T> dataOper;
                    dataOper.setAcquiredPoints(acq);

                    gtPlusWaveletNoNullSpace2DOperator<T> wavNoNullSpace2DOperator;
                    wavNoNullSpace2DOperator.setAcquiredPoints(acq);

                    if ( workOrder2DT->spirit_use_coil_sen_map_ && workOrder2DT->coilMap_ )
                    {
                        wavNoNullSpace2DOperator.setCoilSenMap(coilMapS);
                    }

                    ncgsolver.add(spirit_noNullSpace, T( (value_type)workOrder2DT->spirit_parallel_imaging_lamda_ ) );
                    ncgsolver.add(wavNoNullSpace2DOperator, T( (value_type)workOrder2DT->spirit_image_reg_lamda_ ) );
                    ncgsolver.add(dataOper, T( (value_type)workOrder2DT->spirit_data_fidelity_lamda_ ) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 2D without null space ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_2D_res_noNullSpace"); }
                }
            }

            Gadgetron::scal(T( (value_type)workOrder2DT->spirit_ncg_scale_factor_ ), res);
        }
        else
        {
            res = kspaceLinear;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DTL1SPIRITNCG<T>::performUnwarppingImpl(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DTL1SPIRITNCG<T>::
performUnwarppingImpl(gtPlusReconJob2DT<T>& job)
{
    try
    {
        hoNDArray<T>& kspace = job.kspace;
        hoNDArray<T>& ker = job.ker;
        hoNDArray<T>& res = job.res;
        gtPlusReconWorkOrder<T>* workOrder2DT = &(job.workOrder2DT);

        GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(workOrder2DT, kspace, ker, res, job.job_index_S_));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DTL1SPIRITNCG<T>::performUnwarppingImpl(job) ... ");
        return false;
    }

    return true;
}

}}
