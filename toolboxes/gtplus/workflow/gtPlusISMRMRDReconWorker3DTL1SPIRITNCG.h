/** \file   gtPlusISMRMRDReconWorker3DTL1SPIRITNCG.h
    \brief  Implement the 3DT non-linear SPIRIT reconstruction using the non-linear CG solver
    \author Hui Xue
*/

#pragma once

#include "gtPlusISMRMRDReconWorker3DTSPIRIT.h"
#include "gtPlusSPIRIT2DTOperator.h"
#include "gtPlusSPIRITNoNullSpace2DOperator.h"
#include "gtPlusSPIRITNoNullSpace2DTOperator.h"
#include "gtPlusNCGSolver.h"
#include "gtPlusWavelet2DOperator.h"
#include "gtPlusWavelet3DOperator.h"
#include "gtPlusWaveletNoNullSpace2DOperator.h"
#include "gtPlusWaveletNoNullSpace3DOperator.h"
#include "gtPlusDataFidelityOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorker3DTL1SPIRITNCG : public gtPlusReconWorker3DTSPIRIT<T>
{
public:

    typedef gtPlusReconWorker3DTSPIRIT<T> BaseClass;
    typedef gtPlusReconWorkOrder3DT<T> WorkOrderType;
    typedef typename BaseClass::value_type value_type;

    gtPlusReconWorker3DTL1SPIRITNCG() : BaseClass() {}
    virtual ~gtPlusReconWorker3DTL1SPIRITNCG() {}

    virtual bool performUnwarppingImpl(gtPlusReconWorkOrder<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res, size_t n);
    virtual bool performUnwarppingImplROPermuted(gtPlusReconWorkOrder<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& kernel, hoNDArray<T>& coilMap, hoNDArray<T>& res);
    virtual bool performUnwarppingImplROPermuted(gtPlusReconWorkOrder<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& kernel, hoNDArray<T>& coilMap, hoNDArray<T>& kspaceLinear, hoNDArray<T>& res);
    virtual bool performUnwarppingImpl(gtPlusReconJob2DT<T>& job);

    virtual bool computeKSpace(gtPlusReconWorkOrder3DT<T>* workOrder3DT);

    virtual bool autoReconParameter(gtPlusReconWorkOrder<T>* workOrder);

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::verbose_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_cplx_;
    using BaseClass::gtPlus_mem_manager_;

//protected::

    using BaseClass::ref_src_;
    using BaseClass::ref_dst_;
    using BaseClass::data_dst_;
    using BaseClass::ref_coil_map_dst_;
    using BaseClass::startE1_;
    using BaseClass::endE1_;

    gtPlusSPIRIT<T> spirit_;
};

template <typename T> 
bool gtPlusReconWorker3DTL1SPIRITNCG<T>::computeKSpace(gtPlusReconWorkOrder3DT<T>* workOrder3DT)
{
    bool recon_kspace = true;
    if ( workOrder3DT->spirit_perform_nonlinear_ && workOrder3DT->spirit_use_coil_sen_map_ ) recon_kspace = false;
    return recon_kspace;
}

template <typename T> 
bool gtPlusReconWorker3DTL1SPIRITNCG<T>::autoReconParameter(gtPlusReconWorkOrder<T>* workOrder)
{
    BaseClass::autoReconParameter(workOrder);

    gtPlusReconWorkOrder3DT<T>* workOrder3DT = dynamic_cast<gtPlusReconWorkOrder3DT<T>*>(workOrder);
    if ( workOrder3DT == NULL ) return false;

    double acceFactor = workOrder3DT->acceFactorE1_ * workOrder3DT->acceFactorE2_;

    if ( workOrder3DT->spirit_perform_linear_ )
    {
        if ( acceFactor>=16 )
        {
            workOrder3DT->spirit_3D_scale_per_chunk_ = true;

            if ( workOrder3DT->spirit_solve_symmetric_ )
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.0025;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
            else
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.0025;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
        }
        else if ( acceFactor>=12 )
        {
            workOrder3DT->spirit_3D_scale_per_chunk_ = true;

            if ( workOrder3DT->spirit_solve_symmetric_ )
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.0025;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
            else
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.0025;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
        }
        else if ( acceFactor>=9 )
        {
            workOrder3DT->spirit_3D_scale_per_chunk_ = true;

            if ( workOrder3DT->spirit_solve_symmetric_ )
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.0025;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
            else
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.0025;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
        }
        else if ( acceFactor>=6 )
        {
            workOrder3DT->spirit_3D_scale_per_chunk_ = true;

            if ( workOrder3DT->spirit_solve_symmetric_ )
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.002;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
            else
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.002;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
        }
        else if ( acceFactor>=4 )
        {
            workOrder3DT->spirit_3D_scale_per_chunk_ = true;

            if ( workOrder3DT->spirit_solve_symmetric_ )
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.0015;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
            else
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.002;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
        }
        else
        {
            workOrder3DT->spirit_3D_scale_per_chunk_ = true;

            if ( workOrder3DT->spirit_solve_symmetric_ )
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.0015;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
            else
            {
                workOrder3DT->spirit_image_reg_lamda_ = 0.002;
                workOrder3DT->spirit_ncg_iter_thres_ = 0.001;
            }
        }
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTL1SPIRITNCG<T>::
performUnwarppingImpl(gtPlusReconWorkOrder<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res, size_t n)
{
    try
    {
        // RO, E1, E2, srcCHA, dstCHA

        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t E2 = kspace.get_size(2);

        size_t srcCHA = adj_forward_G_I.get_size(3);
        size_t dstCHA = adj_forward_G_I.get_size(4);

        res.create(kspace.get_dimensions());

        // perform the 3D recon by read-out decoupling

        hoNDArrayMemoryManaged<T> kspaceIfftRO(RO, E1, E2, srcCHA, gtPlus_mem_manager_);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(kspace, kspaceIfftRO);
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIfftRO, debugFolder_+"kspaceIfftRO"); }

        hoNDArrayMemoryManaged<T> kspaceIfftROPermuted(E1, E2, srcCHA, RO, gtPlus_mem_manager_);

        if ( performTiming_ ) { gt_timer3_.start("permtue RO to 4th dimension ... "); }
        GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo4thDimensionFor3DRecon(kspaceIfftRO, kspaceIfftROPermuted));
        if ( performTiming_ ) { gt_timer3_.stop(); }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIfftROPermuted, debugFolder_+"kspaceIfftROPermuted"); }

        // permute kernel
        hoNDArray<T> kerPermuted(E1, E2, srcCHA, dstCHA, RO);
        if ( performTiming_ ) { gt_timer3_.start("permute kernel RO to 5th dimension ... "); }
        GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteE2To5thDimension( adj_forward_G_I, kerPermuted));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        // permute coil map
        hoNDArray<T> coilMapN(RO, E1, E2, dstCHA, workOrder3DT->coilMap_->begin()+n*RO*E1*E2*dstCHA);
        hoNDArray<T> coilMapPermuted(E1, E2, dstCHA, RO);
        if ( performTiming_ ) { gt_timer3_.start("permtue coil map RO to 4th dimension ... "); }
        GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo4thDimensionFor3DRecon(coilMapN, coilMapPermuted));
        if ( performTiming_ ) { gt_timer3_.stop(); }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(coilMapPermuted, debugFolder_+"coilMapPermuted"); }

        hoNDArray<T> resPermuted(E1, E2, dstCHA, RO);
        GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImplROPermuted(workOrder3DT, kspaceIfftROPermuted, kerPermuted, coilMapPermuted, resPermuted));

        // permute the unwrapped kspace
        if ( performTiming_ ) { gt_timer3_.start("permtue RO to 1st dimension ... "); }
        GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo1stDimensionFor3DRecon(resPermuted, kspaceIfftRO));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        // perform fft along the first dimension
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(kspaceIfftRO, res);
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"res_3DSpirit"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTL1SPIRITNCG<T>::performUnwarppingImpl(gtPlusReconWorkOrder3DT<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTL1SPIRITNCG<T>::
performUnwarppingImplROPermuted(gtPlusReconWorkOrder<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& kernel, hoNDArray<T>& coilMap, hoNDArray<T>& res)
{
    try
    {
        size_t E1 = kspace.get_size(0);
        size_t E2 = kspace.get_size(1);
        size_t RO = kspace.get_size(3);

        size_t kerE1 = kernel.get_size(0);
        size_t kerE2 = kernel.get_size(1);
        size_t srcCHA = kernel.get_size(2);
        size_t dstCHA = kernel.get_size(3);
        size_t kerN = kernel.get_size(5);

        hoNDArray<T>* kerIm = &kernel;
        hoNDArray<T> kerImE1E2RO;
        if ( kerE1!=E1 || kerE2!=E2 )
        {
            GDEBUG_STREAM("gtPlusReconWorker3DTL1SPIRITNCG, kerE1!=E1 || kerE2!=E2, kernel needs to be converted along E1 and E2 ... ");

            if ( gtPlus_mem_manager_ )
            {
                // kerImE1E2RO will be cleared as all '0' 
                kerImE1E2RO.create(E1, E2, srcCHA, dstCHA, RO, kerN, (T*)(gtPlus_mem_manager_->allocate(sizeof(T)*(size_t)RO*E1*E2*srcCHA*dstCHA)));
            }
            else
            {
                kerImE1E2RO.create(E1, E2, srcCHA, dstCHA, RO, kerN);
                Gadgetron::clear(kerImE1E2RO);
            }

            GADGET_CHECK_RETURN_FALSE(spirit_.imageDomainKernelE1E2RO(kernel, (int)E1, (int)E2, kerImE1E2RO));
            kerIm = &kerImE1E2RO;
        }

        hoNDArray<T> kspaceLinear(kspace);
        res = kspace;

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace"); }

        bool performLinear = workOrder3DT->spirit_perform_linear_;
        if ( !workOrder3DT->spirit_perform_nonlinear_ ) performLinear = true;

        if ( performLinear )
        {
            if ( performTiming_ ) { gt_timer3_.start("NCG spirit linear solver for 3DT ... "); }
            GADGET_CHECK_RETURN_FALSE(BaseClass::performUnwarppingImplROPermuted(workOrder3DT, kspace, *kerIm, coilMap, kspaceLinear));
            if ( performTiming_ ) { gt_timer3_.stop(); }
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceLinear, debugFolder_+"kspaceLinear"); }

        if ( workOrder3DT->spirit_perform_nonlinear_ )
        {
            if ( workOrder3DT->spirit_3D_scale_per_chunk_ )
            {
                typename realType<T>::Type scaleFactor = 1.0;
                Gadgetron::norm2(kspace, scaleFactor);
                scaleFactor /= (value_type)(RO*std::sqrt(double(srcCHA)));

                workOrder3DT->spirit_ncg_scale_factor_ = scaleFactor;
            }

            // apply the scale
            Gadgetron::scal( (value_type)(1.0/workOrder3DT->spirit_ncg_scale_factor_), kspaceLinear);
            Gadgetron::scal( (value_type)(1.0/workOrder3DT->spirit_ncg_scale_factor_), kspace);

            boost::shared_ptr< hoNDArray<T> > coilMapN;
            if ( workOrder3DT->coilMap_ 
                && workOrder3DT->coilMap_->get_size(0)==E1 
                && workOrder3DT->coilMap_->get_size(1)==E2 
                && workOrder3DT->coilMap_->get_size(2)==dstCHA 
                && workOrder3DT->coilMap_->get_size(3)==RO )
            {
                coilMapN = boost::shared_ptr< hoNDArray<T> >( new hoNDArray<T>(E1, E2, dstCHA, RO, coilMap.begin()) );
            }

            if ( RO > 1 )
            {
                boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(E1, E2, srcCHA, dstCHA, RO, kerIm->begin()));
                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(E1, E2, srcCHA, RO, kspace.begin()));

                gtPlusNCGSolver<hoNDArray<T>, hoNDArray<T>, gtPlusOperator<T> > ncgsolver;
                ncgsolver.iterMax_ = workOrder3DT->spirit_ncg_iter_max_;
                ncgsolver.printIter_ = workOrder3DT->spirit_ncg_print_iter_;
                ncgsolver.secantRatio_ = 2;
                ncgsolver.x0_ = &kspaceLinear;

                hoNDArray<T> b;

                if ( workOrder3DT->spirit_data_fidelity_lamda_ <= 0 )
                {
                    // parallel imaging term
                    gtPlusSPIRIT2DTOperator<T> spirit;
                    spirit.use_symmetric_spirit_ = false;
                    spirit.setMemoryManager(gtPlus_mem_manager_);

                    spirit.setForwardKernel(ker, true);
                    spirit.setAcquiredPoints(acq);

                    // L1 term
                    gtPlusWavelet3DOperator<T> wavNullSpace3DOperator;
                    wavNullSpace3DOperator.setMemoryManager(gtPlus_mem_manager_);
                    wavNullSpace3DOperator.setAcquiredPoints(acq);

                    wavNullSpace3DOperator.scale_factor_first_dimension_ = (value_type)workOrder3DT->spirit_E1_enhancement_ratio_;
                    wavNullSpace3DOperator.scale_factor_second_dimension_ = (value_type)workOrder3DT->spirit_E2_enhancement_ratio_;
                    wavNullSpace3DOperator.scale_factor_third_dimension_ = (value_type)workOrder3DT->spirit_RO_enhancement_ratio_;

                    if ( workOrder3DT->spirit_use_coil_sen_map_ && coilMapN )
                    {
                        wavNullSpace3DOperator.setCoilSenMap(coilMapN);
                    }

                    // set operators
                    ncgsolver.add(spirit, (value_type)(workOrder3DT->spirit_parallel_imaging_lamda_) );
                    ncgsolver.add(wavNullSpace3DOperator, (value_type)(workOrder3DT->spirit_image_reg_lamda_) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 3DT ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_3DT_res"); }

                    spirit.restoreAcquiredKSpace(kspace, res);

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_3DT_res_restored"); }
                }
                else
                {
                    gtPlusSPIRITNoNullSpace2DTOperator<T> spirit_noNullSpace;
                    spirit_noNullSpace.use_symmetric_spirit_ = false;
                    spirit_noNullSpace.setMemoryManager(gtPlus_mem_manager_);

                    spirit_noNullSpace.setForwardKernel(ker, true);
                    spirit_noNullSpace.setAcquiredPoints(acq);

                    gtPlusDataFidelityOperator<T> dataOper;
                    dataOper.setAcquiredPoints(acq);

                    gtPlusWaveletNoNullSpace3DOperator<T> wavNoNullSpace3DOperator;
                    wavNoNullSpace3DOperator.setMemoryManager(gtPlus_mem_manager_);
                    wavNoNullSpace3DOperator.setAcquiredPoints(acq);

                    wavNoNullSpace3DOperator.scale_factor_first_dimension_ = (value_type)workOrder3DT->spirit_E1_enhancement_ratio_;
                    wavNoNullSpace3DOperator.scale_factor_second_dimension_ = (value_type)workOrder3DT->spirit_E2_enhancement_ratio_;
                    wavNoNullSpace3DOperator.scale_factor_third_dimension_ = (value_type)workOrder3DT->spirit_RO_enhancement_ratio_;

                    if ( workOrder3DT->spirit_use_coil_sen_map_ && coilMapN )
                    {
                        wavNoNullSpace3DOperator.setCoilSenMap(coilMapN);
                    }

                    ncgsolver.add(spirit_noNullSpace, (value_type)(workOrder3DT->spirit_parallel_imaging_lamda_) );
                    ncgsolver.add(wavNoNullSpace3DOperator, (value_type)(workOrder3DT->spirit_image_reg_lamda_) );
                    ncgsolver.add(dataOper, (value_type)(workOrder3DT->spirit_data_fidelity_lamda_) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 3DT without null space ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_3DT_res_noNullSpace"); }
                }
            }
            else
            {
                boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(E1, E2, srcCHA, dstCHA, kerIm->begin()));
                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(E1, E2, srcCHA, kspace.begin()));

                gtPlusNCGSolver<hoNDArray<T>, hoNDArray<T>, gtPlusOperator<T> > ncgsolver;
                ncgsolver.iterMax_ = workOrder3DT->spirit_ncg_iter_max_;
                ncgsolver.printIter_ = workOrder3DT->spirit_ncg_print_iter_;
                ncgsolver.secantRatio_ = 2;
                ncgsolver.x0_ = &kspaceLinear;

                hoNDArray<T> b;

                if ( workOrder3DT->spirit_data_fidelity_lamda_ <= 0 )
                {
                    // parallel imaging term
                    gtPlusSPIRIT2DOperator<T> spirit;
                    spirit.use_symmetric_spirit_ = false;
                    spirit.setMemoryManager(gtPlus_mem_manager_);
                    spirit.setForwardKernel(ker, true);
                    spirit.setAcquiredPoints(acq);

                    // L1 term
                    gtPlusWavelet2DOperator<T> wavNullSpace2DOperator;
                    wavNullSpace2DOperator.setMemoryManager(gtPlus_mem_manager_);
                    wavNullSpace2DOperator.setAcquiredPoints(acq);

                    if ( workOrder3DT->spirit_use_coil_sen_map_ && coilMapN )
                    {
                        wavNullSpace2DOperator.setCoilSenMap(coilMapN);
                    }

                    // set operators
                    ncgsolver.add(spirit, (value_type)(workOrder3DT->spirit_parallel_imaging_lamda_) );
                    ncgsolver.add(wavNullSpace2DOperator, (value_type)(workOrder3DT->spirit_image_reg_lamda_) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 3D ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_3D_res"); }

                    spirit.restoreAcquiredKSpace(kspace, res);

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_3D_res_restored"); }
                }
                else
                {
                    gtPlusSPIRITNoNullSpace2DOperator<T> spirit_noNullSpace;
                    spirit_noNullSpace.use_symmetric_spirit_ = false;
                    spirit_noNullSpace.setMemoryManager(gtPlus_mem_manager_);
                    spirit_noNullSpace.setForwardKernel(ker, true);
                    spirit_noNullSpace.setAcquiredPoints(acq);

                    gtPlusDataFidelityOperator<T> dataOper;
                    dataOper.setAcquiredPoints(acq);

                    gtPlusWaveletNoNullSpace2DOperator<T> wavNoNullSpace2DOperator;
                    wavNoNullSpace2DOperator.setMemoryManager(gtPlus_mem_manager_);
                    wavNoNullSpace2DOperator.setAcquiredPoints(acq);

                    if ( workOrder3DT->spirit_use_coil_sen_map_ && coilMapN )
                    {
                        wavNoNullSpace2DOperator.setCoilSenMap(coilMapN);
                    }

                    ncgsolver.add(spirit_noNullSpace, (value_type)(workOrder3DT->spirit_parallel_imaging_lamda_) );
                    ncgsolver.add(wavNoNullSpace2DOperator, (value_type)(workOrder3DT->spirit_image_reg_lamda_) );
                    ncgsolver.add(dataOper, (value_type)(workOrder3DT->spirit_data_fidelity_lamda_) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 3D without null space ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"ncg_spirit_3D_res_noNullSpace"); }
                }
            }

            Gadgetron::scal( (value_type)(workOrder3DT->spirit_ncg_scale_factor_), res);
        }
        else
        {
            res = kspaceLinear;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTL1SPIRITNCG<T>::performUnwarppingImplROPermuted(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTL1SPIRITNCG<T>::
performUnwarppingImplROPermuted(gtPlusReconWorkOrder<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& kernel, hoNDArray<T>& coilMap, hoNDArray<T>& kspaceLinear, hoNDArray<T>& res)
{
    try
    {
        size_t E1 = kspace.get_size(0);
        size_t E2 = kspace.get_size(1);
        size_t RO = kspace.get_size(3);

        size_t kerE1 = kernel.get_size(0);
        size_t kerE2 = kernel.get_size(1);
        size_t srcCHA = kernel.get_size(2);
        size_t dstCHA = kernel.get_size(3);
        size_t kerN = kernel.get_size(5);

        hoNDArray<T>* kerIm = &kernel;
        hoNDArray<T> kerImE1E2RO;
        if ( kerE1!=E1 || kerE2!=E2 )
        {
            GDEBUG_STREAM("gtPlusReconWorker3DTL1SPIRITNCG, kerE1!=E1 || kerE2!=E2, kernel needs to be converted along E1 and E2 ... ");

            if ( gtPlus_mem_manager_ )
            {
                // kerImE1E2RO will be cleared as all '0' 
                kerImE1E2RO.create(E1, E2, srcCHA, dstCHA, RO, kerN, (T*)(gtPlus_mem_manager_->allocate(sizeof(T)*(size_t)RO*E1*E2*srcCHA*dstCHA)));
            }
            else
            {
                kerImE1E2RO.create(E1, E2, srcCHA, dstCHA, RO, kerN);
                Gadgetron::clear(kerImE1E2RO);
            }

            GADGET_CHECK_RETURN_FALSE(spirit_.imageDomainKernelE1E2RO(kernel, (int)E1, (int)E2, kerImE1E2RO));
            kerIm = &kerImE1E2RO;
        }

        res = kspace;

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace"); }

        bool performLinear = workOrder3DT->spirit_perform_linear_;
        if ( !workOrder3DT->spirit_perform_nonlinear_ ) performLinear = true;

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceLinear, debugFolder_+"kspaceLinear"); }

        if ( workOrder3DT->spirit_perform_nonlinear_ )
        {
            if ( workOrder3DT->spirit_3D_scale_per_chunk_ )
            {
                typename realType<T>::Type scaleFactor = 1.0;
                Gadgetron::norm2(kspace, scaleFactor);
                scaleFactor /= (value_type)(RO*std::sqrt(double(srcCHA)));

                workOrder3DT->spirit_ncg_scale_factor_ = scaleFactor;
            }

            // apply the scale
            Gadgetron::scal((value_type)(1.0/workOrder3DT->spirit_ncg_scale_factor_), kspaceLinear);
            Gadgetron::scal((value_type)(1.0/workOrder3DT->spirit_ncg_scale_factor_), kspace);

            boost::shared_ptr< hoNDArray<T> > coilMapN;
            if ( workOrder3DT->coilMap_ 
                && workOrder3DT->coilMap_->get_size(0)==E1 
                && workOrder3DT->coilMap_->get_size(1)==E2 
                && workOrder3DT->coilMap_->get_size(2)==dstCHA 
                && workOrder3DT->coilMap_->get_size(3)==RO )
            {
                coilMapN = boost::shared_ptr< hoNDArray<T> >( new hoNDArray<T>(E1, E2, dstCHA, RO, coilMap.begin()) );
            }

            if ( RO > 1 )
            {
                boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(E1, E2, srcCHA, dstCHA, RO, kerIm->begin()));
                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(E1, E2, srcCHA, RO, kspace.begin()));

                gtPlusNCGSolver<hoNDArray<T>, hoNDArray<T>, gtPlusOperator<T> > ncgsolver;
                ncgsolver.iterMax_ = workOrder3DT->spirit_ncg_iter_max_;
                ncgsolver.printIter_ = workOrder3DT->spirit_ncg_print_iter_;
                ncgsolver.secantRatio_ = 2;
                ncgsolver.x0_ = &kspaceLinear;

                hoNDArray<T> b;

                if ( workOrder3DT->spirit_data_fidelity_lamda_ <= 0 )
                {
                    // parallel imaging term
                    gtPlusSPIRIT2DTOperator<T> spirit;
                    spirit.use_symmetric_spirit_ = false;
                    spirit.setMemoryManager(gtPlus_mem_manager_);

                    spirit.setForwardKernel(ker, true);
                    spirit.setAcquiredPoints(acq);

                    // L1 term
                    gtPlusWavelet3DOperator<T> wavNullSpace3DOperator;
                    wavNullSpace3DOperator.setMemoryManager(gtPlus_mem_manager_);
                    wavNullSpace3DOperator.setAcquiredPoints(acq);

                    wavNullSpace3DOperator.scale_factor_first_dimension_ = (value_type)workOrder3DT->spirit_E1_enhancement_ratio_;
                    wavNullSpace3DOperator.scale_factor_second_dimension_ = (value_type)workOrder3DT->spirit_E2_enhancement_ratio_;
                    wavNullSpace3DOperator.scale_factor_third_dimension_ = (value_type)workOrder3DT->spirit_RO_enhancement_ratio_;

                    if ( workOrder3DT->spirit_use_coil_sen_map_ && coilMapN )
                    {
                        wavNullSpace3DOperator.setCoilSenMap(coilMapN);
                    }

                    // set operators
                    ncgsolver.add(spirit, (value_type)(workOrder3DT->spirit_parallel_imaging_lamda_) );
                    ncgsolver.add(wavNullSpace3DOperator, (value_type)(workOrder3DT->spirit_image_reg_lamda_) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 3DT ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, "ncg_spirit_3DT_res"); }

                    spirit.restoreAcquiredKSpace(kspace, res);

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, "ncg_spirit_3DT_res_restored"); }
                }
                else
                {
                    gtPlusSPIRITNoNullSpace2DTOperator<T> spirit_noNullSpace;
                    spirit_noNullSpace.use_symmetric_spirit_ = false;
                    spirit_noNullSpace.setMemoryManager(gtPlus_mem_manager_);

                    spirit_noNullSpace.setForwardKernel(ker, true);
                    spirit_noNullSpace.setAcquiredPoints(acq);

                    gtPlusDataFidelityOperator<T> dataOper;
                    dataOper.setAcquiredPoints(acq);

                    gtPlusWaveletNoNullSpace3DOperator<T> wavNoNullSpace3DOperator;
                    wavNoNullSpace3DOperator.setMemoryManager(gtPlus_mem_manager_);
                    wavNoNullSpace3DOperator.setAcquiredPoints(acq);

                    wavNoNullSpace3DOperator.scale_factor_first_dimension_ = (value_type)workOrder3DT->spirit_E1_enhancement_ratio_;
                    wavNoNullSpace3DOperator.scale_factor_second_dimension_ = (value_type)workOrder3DT->spirit_E2_enhancement_ratio_;
                    wavNoNullSpace3DOperator.scale_factor_third_dimension_ = (value_type)workOrder3DT->spirit_RO_enhancement_ratio_;

                    if ( workOrder3DT->spirit_use_coil_sen_map_ && coilMapN )
                    {
                        wavNoNullSpace3DOperator.setCoilSenMap(coilMapN);
                    }

                    ncgsolver.add(spirit_noNullSpace, (value_type)(workOrder3DT->spirit_parallel_imaging_lamda_) );
                    ncgsolver.add(wavNoNullSpace3DOperator, (value_type)(workOrder3DT->spirit_image_reg_lamda_) );
                    ncgsolver.add(dataOper, (value_type)(workOrder3DT->spirit_data_fidelity_lamda_) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 3DT without null space ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, "ncg_spirit_3DT_res_noNullSpace"); }
                }
            }
            else
            {
                boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(E1, E2, srcCHA, dstCHA, kerIm->begin()));
                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(E1, E2, srcCHA, kspace.begin()));

                gtPlusNCGSolver<hoNDArray<T>, hoNDArray<T>, gtPlusOperator<T> > ncgsolver;
                ncgsolver.iterMax_ = workOrder3DT->spirit_ncg_iter_max_;
                ncgsolver.printIter_ = workOrder3DT->spirit_ncg_print_iter_;
                ncgsolver.secantRatio_ = 2;
                ncgsolver.x0_ = &kspaceLinear;

                hoNDArray<T> b;

                if ( workOrder3DT->spirit_data_fidelity_lamda_ <= 0 )
                {
                    // parallel imaging term
                    gtPlusSPIRIT2DOperator<T> spirit;
                    spirit.use_symmetric_spirit_ = false;
                    spirit.setMemoryManager(gtPlus_mem_manager_);
                    spirit.setForwardKernel(ker, true);
                    spirit.setAcquiredPoints(acq);

                    // L1 term
                    gtPlusWavelet2DOperator<T> wavNullSpace2DOperator;
                    wavNullSpace2DOperator.setMemoryManager(gtPlus_mem_manager_);
                    wavNullSpace2DOperator.setAcquiredPoints(acq);

                    if ( workOrder3DT->spirit_use_coil_sen_map_ && coilMapN )
                    {
                        wavNullSpace2DOperator.setCoilSenMap(coilMapN);
                    }

                    // set operators
                    ncgsolver.add(spirit, (value_type)(workOrder3DT->spirit_parallel_imaging_lamda_) );
                    ncgsolver.add(wavNullSpace2DOperator, (value_type)(workOrder3DT->spirit_image_reg_lamda_) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 3D ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, "ncg_spirit_3D_res"); }

                    spirit.restoreAcquiredKSpace(kspace, res);

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, "ncg_spirit_3D_res_restored"); }
                }
                else
                {
                    gtPlusSPIRITNoNullSpace2DOperator<T> spirit_noNullSpace;
                    spirit_noNullSpace.use_symmetric_spirit_ = false;
                    spirit_noNullSpace.setMemoryManager(gtPlus_mem_manager_);
                    spirit_noNullSpace.setForwardKernel(ker, true);
                    spirit_noNullSpace.setAcquiredPoints(acq);

                    gtPlusDataFidelityOperator<T> dataOper;
                    dataOper.setAcquiredPoints(acq);

                    gtPlusWaveletNoNullSpace2DOperator<T> wavNoNullSpace2DOperator;
                    wavNoNullSpace2DOperator.setMemoryManager(gtPlus_mem_manager_);
                    wavNoNullSpace2DOperator.setAcquiredPoints(acq);

                    if ( workOrder3DT->spirit_use_coil_sen_map_ && coilMapN )
                    {
                        wavNoNullSpace2DOperator.setCoilSenMap(coilMapN);
                    }

                    ncgsolver.add(spirit_noNullSpace, (value_type)(workOrder3DT->spirit_parallel_imaging_lamda_) );
                    ncgsolver.add(wavNoNullSpace2DOperator, (value_type)(workOrder3DT->spirit_image_reg_lamda_) );
                    ncgsolver.add(dataOper, (value_type)(workOrder3DT->spirit_data_fidelity_lamda_) );

                    if ( performTiming_ ) { gt_timer3_.start("NCG spirit solver for 3D without null space ... "); }
                    ncgsolver.solve(b, res);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, "ncg_spirit_3D_res_noNullSpace"); }
                }
            }

            Gadgetron::scal( (value_type)(workOrder3DT->spirit_ncg_scale_factor_), res);
        }
        else
        {
            res = kspaceLinear;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTL1SPIRITNCG<T>::performUnwarppingImplROPermuted(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTL1SPIRITNCG<T>::
performUnwarppingImpl(gtPlusReconJob2DT<T>& job)
{
    try
    {
        hoNDArray<T>& kspace = job.kspace; // [E1 E2 srcCHA RO 1]
        hoNDArray<T>& ker = job.ker; // [E1 E2 srcCHA dstCHA RO 1]
        hoNDArray<T>& res = job.res; // [E1 E2 dstCHA RO 1]
        gtPlusReconWorkOrder<T>* workOrder3DT = &(job.workOrder2DT);

        job.res = job.kspace;

        GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImplROPermuted(workOrder3DT, kspace, ker, *job.workOrder2DT.coilMap_, res));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTL1SPIRITNCG<T>::performUnwarppingImpl(job) ... ");
        return false;
    }

    return true;
}

}}
