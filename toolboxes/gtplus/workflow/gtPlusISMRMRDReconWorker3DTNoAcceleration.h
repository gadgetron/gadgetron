/** \file   gtPlusISMRMRDReconWorker3DTNoAcceleration.h
    \brief  Implement the 3DT reconstruction without the k-space undersampling
    \author Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"

#include "GadgetronTimer.h"

#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorker3DT.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorker3DTNoAcceleration : public gtPlusReconWorker3DT<T>
{
public:

    typedef gtPlusReconWorker3DT<T> BaseClass;
    typedef typename BaseClass::value_type value_type;

    gtPlusReconWorker3DTNoAcceleration() : BaseClass() {}
    virtual ~gtPlusReconWorker3DTNoAcceleration() {}

    virtual bool performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT);

    virtual bool computeKSpace(gtPlusReconWorkOrder3DT<T>* /*workOrder3DT*/) { return false; }

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::verbose_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_cplx_;

    using BaseClass::ref_src_;
    using BaseClass::ref_dst_;
    using BaseClass::data_dst_;
    using BaseClass::ref_coil_map_dst_;
    using BaseClass::startE1_;
    using BaseClass::endE1_;
    using BaseClass::startE2_;
    using BaseClass::endE2_;
};

template <typename T> 
bool gtPlusReconWorker3DTNoAcceleration<T>::performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(workOrder3DT!=NULL);

        if ( !workOrder3DT->workFlow_use_BufferedKernel_ )
        {
            if ( performTiming_ ) { gt_timer1_.start("prepRef"); }
            GADGET_CHECK_RETURN_FALSE(this->prepRef(workOrder3DT, workOrder3DT->ref_, 
                                                workOrder3DT->ref_recon_, 
                                                workOrder3DT->ref_coil_map_, 
                                                workOrder3DT->start_RO_, workOrder3DT->end_RO_, 
                                                workOrder3DT->start_E1_, workOrder3DT->end_E1_, 
                                                workOrder3DT->start_E2_, workOrder3DT->end_E2_, 
                                                workOrder3DT->data_.get_size(1), workOrder3DT->data_.get_size(2)));
            if ( performTiming_ ) { gt_timer1_.stop(); }
        }

        size_t RO = workOrder3DT->data_.get_size(0);
        size_t E1 = workOrder3DT->data_.get_size(1);
        size_t E2 = workOrder3DT->data_.get_size(2);
        size_t CHA = workOrder3DT->data_.get_size(3);
        size_t N = workOrder3DT->data_.get_size(4);

        size_t refN = workOrder3DT->ref_recon_.get_size(4);
        size_t usedN;

        // estimate the coil sensitivity
        if ( !workOrder3DT->workFlow_use_BufferedKernel_ 
                    || (workOrder3DT->coilMap_->get_size(0)!=RO) 
                    || (workOrder3DT->coilMap_->get_size(1)!=E1)
                    || (workOrder3DT->coilMap_->get_size(2)!=E2) )
        {
            workOrder3DT->coilMap_->create(RO, E1, E2, CHA, refN);

            if ( workOrder3DT->no_acceleration_same_combinationcoeff_allN_ )
            {
                usedN = workOrder3DT->no_acceleration_whichN_combinationcoeff_;
                if ( usedN >= refN ) usedN = refN-1;

                hoNDArray<T> refCoilMapN(RO, E1, E2, CHA, workOrder3DT->ref_coil_map_.begin()+usedN*RO*E1*E2*CHA);

                hoNDArray<T> buffer3DT(refCoilMapN.get_dimensions());

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(refCoilMapN, buffer3DT);

                hoNDArray<T> coilMapN(RO, E1, E2, CHA, workOrder3DT->coilMap_->begin()+usedN*RO*E1*E2*CHA);

                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(buffer3DT, 
                        coilMapN, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));

                GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder3DT->coilMap_, usedN));
            }
            else
            {
                hoNDArray<T> buffer3DT(workOrder3DT->ref_coil_map_.get_dimensions());

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->ref_coil_map_, buffer3DT);

                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(buffer3DT, 
                        *workOrder3DT->coilMap_, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder3DT->coilMap_, debugFolder_+"coilMap_"); }
        }

        size_t num = 0;

        size_t e1, e2, n;
        for (n = 0; n < N; n++)
        {
            for (e2 = 0; e2 < E2; e2++)
            {
                for (e1 = 0; e1 < E1; e1++)
                {
                    if (std::abs(workOrder3DT->data_(RO / 2, e1, e2, 0, n)) > 0)
                    {
                        num++;
                    }
                }
            }
        }

        if (num > 0)
        {
            double lenRO = RO;
            if (!(workOrder3DT->start_RO_<0 || workOrder3DT->end_RO_<0 || (workOrder3DT->end_RO_ - workOrder3DT->start_RO_ + 1 == RO)))
            {
                lenRO = (workOrder3DT->end_RO_ - workOrder3DT->start_RO_ + 1);
            }
            if (this->verbose_) GDEBUG_STREAM("gtPlusReconWorker3DTNoAcceleration, length for RO : " << lenRO);

            double effectiveAcceFactor = (double)(N*E1*E2) / (num);
            if (this->verbose_) GDEBUG_STREAM("gtPlusReconWorker3DTNoAcceleration, effectiveAcceFactor : " << effectiveAcceFactor);

            typename realType<T>::Type fftCompensationRatio = (typename realType<T>::Type)(std::sqrt(RO*effectiveAcceFactor / lenRO));
            if (this->verbose_) GDEBUG_STREAM("gtPlusReconWorker3DTNoAcceleration, fftCompensationRatio : " << fftCompensationRatio);

            Gadgetron::scal(fftCompensationRatio, workOrder3DT->data_);
        }

        // partial fourier handling
        GADGET_CHECK_RETURN_FALSE(this->performPartialFourierHandling(workOrder3DT));

        workOrder3DT->complexIm_.create(RO, E1, E2, N);

        if ( performTiming_ ) { gt_timer1_.start("perform coil combination"); }

        hoNDArray<T> buffer3DT(workOrder3DT->data_.get_dimensions());
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->data_, buffer3DT);
        // gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(buffer3DT, *workOrder3DT->coilMap_, workOrder3DT->complexIm_);
        Gadgetron::coil_combine(buffer3DT, *workOrder3DT->coilMap_, 3, workOrder3DT->complexIm_);

        if ( performTiming_ ) { gt_timer1_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"combined"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTNoAcceleration<T>::performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT) ... ");
        return false;
    }

    return true;
}

}}
