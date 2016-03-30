/** \file   gtPlusISMRMRDReconWorker2DTNoAcceleration.h
    \brief  Implement the 2DT reconstruction without the k-space undersampling
    \author Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"

#include "GadgetronTimer.h"

#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorker2DT.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorker2DTNoAcceleration : public gtPlusReconWorker2DT<T>
{
public:

    typedef gtPlusReconWorker2DT<T> BaseClass;
    typedef typename BaseClass::value_type value_type;

    gtPlusReconWorker2DTNoAcceleration() : BaseClass() {}
    virtual ~gtPlusReconWorker2DTNoAcceleration() {}

    virtual bool performRecon(gtPlusReconWorkOrder2DT<T>* workOrder2DT);

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
};

template <typename T> 
bool gtPlusReconWorker2DTNoAcceleration<T>::performRecon(gtPlusReconWorkOrder2DT<T>* workOrder2DT)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(workOrder2DT!=NULL);

        if ( !workOrder2DT->workFlow_use_BufferedKernel_ )
        {
            if ( performTiming_ ) { gt_timer1_.start("prepRef"); }
            GADGET_CHECK_RETURN_FALSE(this->prepRef(workOrder2DT, workOrder2DT->ref_, workOrder2DT->ref_recon_, workOrder2DT->ref_coil_map_, 
                                                workOrder2DT->start_RO_, workOrder2DT->end_RO_, workOrder2DT->start_E1_, workOrder2DT->end_E1_, workOrder2DT->data_.get_size(1)));
            if ( performTiming_ ) { gt_timer1_.stop(); }
        }

        size_t RO = workOrder2DT->data_.get_size(0);
        size_t E1 = workOrder2DT->data_.get_size(1);
        size_t CHA = workOrder2DT->data_.get_size(2);
        size_t N = workOrder2DT->data_.get_size(3);
        size_t S = workOrder2DT->data_.get_size(4);

        size_t refN = workOrder2DT->ref_recon_.get_size(3);
        size_t usedS;

        // apply coil compression coefficients
        if ( !workOrder2DT->workFlow_use_BufferedKernel_ 
                    || (workOrder2DT->coilMap_->get_size(0)!=RO) 
                    || (workOrder2DT->coilMap_->get_size(1)!=E1)
                    || (workOrder2DT->coilMap_->get_size(4)!=S) )
        {
            workOrder2DT->coilMap_->create(RO, E1, CHA, refN, S);

            // estimate the coil sensitivity
            if ( workOrder2DT->no_acceleration_same_combinationcoeff_allS_ )
            {
                usedS = workOrder2DT->no_acceleration_whichS_combinationcoeff_;
                if ( usedS >= S ) usedS = S-1;

                hoNDArray<T> refCoilMapS(RO, E1, CHA, refN, workOrder2DT->ref_coil_map_.begin()+usedS*RO*E1*CHA*refN);
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(refCoilMapS, buffer2DT_);

                hoNDArray<T> coilMapS(RO, E1, CHA, refN, workOrder2DT->coilMap_->begin()+usedS*RO*E1*CHA*refN);

                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(buffer2DT_, 
                        coilMapS, workOrder2DT->coil_map_algorithm_, workOrder2DT->csm_kSize_, workOrder2DT->csm_powermethod_num_, workOrder2DT->csm_iter_num_, (value_type)workOrder2DT->csm_iter_thres_, workOrder2DT->csm_use_gpu_));

                GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder2DT->coilMap_, usedS));
            }
            else
            {
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(workOrder2DT->ref_coil_map_, buffer2DT_);

                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(buffer2DT_, 
                        *workOrder2DT->coilMap_, workOrder2DT->coil_map_algorithm_, workOrder2DT->csm_kSize_, workOrder2DT->csm_powermethod_num_, workOrder2DT->csm_iter_num_, (value_type)workOrder2DT->csm_iter_thres_, workOrder2DT->csm_use_gpu_));
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder2DT->coilMap_, debugFolder_+"coilMap_"); }
        }

        size_t num = 0;

        size_t e1, n, s;
        for (s = 0; s < S; s++)
        {
            for (n = 0; n < N; n++)
            {
                for (e1 = 0; e1 < E1; e1++)
                {
                    if (std::abs(workOrder2DT->data_(RO / 2, e1, 0, n, s)) > 0)
                    {
                        num++;
                    }
                }
            }
        }

        if (num > 0)
        {
            double lenRO = RO;
            if (!(workOrder2DT->start_RO_<0 || workOrder2DT->end_RO_<0 || (workOrder2DT->end_RO_ - workOrder2DT->start_RO_ + 1 == RO)))
            {
                lenRO = (workOrder2DT->end_RO_ - workOrder2DT->start_RO_ + 1);
            }
            if (this->verbose_) GDEBUG_STREAM("gtPlusReconWorker2DTNoAcceleration, length for RO : " << lenRO);

            double effectiveAcceFactor = (double)(S*N*E1) / (num);
            if (this->verbose_) GDEBUG_STREAM("gtPlusReconWorker2DTNoAcceleration, effectiveAcceFactor : " << effectiveAcceFactor);

            typename realType<T>::Type fftCompensationRatio = (typename realType<T>::Type)(std::sqrt(RO*effectiveAcceFactor/lenRO));
            if (this->verbose_) GDEBUG_STREAM("gtPlusReconWorker2DTNoAcceleration, fftCompensationRatio : " << fftCompensationRatio);

            Gadgetron::scal(fftCompensationRatio, workOrder2DT->data_);
        }

        // partial fourier handling
        GADGET_CHECK_RETURN_FALSE(this->performPartialFourierHandling(workOrder2DT));

        workOrder2DT->complexIm_.create(RO, E1, N, S);

        if ( performTiming_ ) { gt_timer1_.start("perform coil combination"); }

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(workOrder2DT->data_, buffer2DT_unwrapping_);

        if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(buffer2DT_unwrapping_, debugFolder_ + "unCombined"); }

        /*if ( refN == N )
        {*/
            // gtPlusISMRMRDReconUtilComplex<T>().coilCombine(buffer2DT_unwrapping_, *(workOrder2DT->coilMap_), workOrder2DT->complexIm_ );
        Gadgetron::coil_combine(buffer2DT_unwrapping_, *(workOrder2DT->coilMap_), 2, workOrder2DT->complexIm_ );
        /*}
        else
        {
            for ( usedS=0; usedS<S; usedS++ )
            {
                hoNDArray<T> unwarppedIm(RO, E1, CHA, N, buffer2DT_unwrapping_.begin()+usedS*RO*E1*CHA*N);
                hoNDArray<T> combined(RO, E1, N, workOrder2DT->complexIm_.begin()+usedS*RO*E1*N);

                if ( refN == N )
                {
                    hoNDArray<T> coilMap(RO, E1, CHA, refN, workOrder2DT->coilMap_->begin()+usedS*RO*E1*CHA*refN);
                    gtPlusISMRMRDReconUtilComplex<T>().coilCombine(unwarppedIm, coilMap, combined);
                }
                else
                {
                    hoNDArray<T> coilMap(RO, E1, CHA, workOrder2DT->coilMap_->begin()+usedS*RO*E1*CHA*N);
                    gtPlusISMRMRDReconUtilComplex<T>().coilCombine(unwarppedIm, coilMap, combined);
                }
            }
        }*/
        if ( performTiming_ ) { gt_timer1_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->complexIm_, debugFolder_+"combined"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DTNoAcceleration<T>::performRecon(gtPlusReconWorkOrder2DT<T>* workOrder2DT) ... ");
        return false;
    }

    return true;
}

}}
