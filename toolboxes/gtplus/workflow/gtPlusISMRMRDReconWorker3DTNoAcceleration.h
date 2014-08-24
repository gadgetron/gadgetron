/** \file   gtPlusISMRMRDReconWorker3DTNoAcceleration.h
    \brief  Implement the 3DT reconstruction without the k-space undersampling
    \author Hui Xue
*/

#pragma once

#include "ismrmrd.h"

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
    using BaseClass::gtPlus_mem_manager_;

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
            GADGET_CHECK_PERFORM(performTiming_, gt_timer1_.start("prepRef"));
            GADGET_CHECK_RETURN_FALSE(this->prepRef(workOrder3DT, workOrder3DT->ref_, 
                                                workOrder3DT->ref_recon_, 
                                                workOrder3DT->ref_coil_map_, 
                                                workOrder3DT->start_RO_, workOrder3DT->end_RO_, 
                                                workOrder3DT->start_E1_, workOrder3DT->end_E1_, 
                                                workOrder3DT->start_E2_, workOrder3DT->end_E2_, 
                                                workOrder3DT->data_.get_size(1), workOrder3DT->data_.get_size(2)));
            GADGET_CHECK_PERFORM(performTiming_, gt_timer1_.stop());
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

                hoNDArrayMemoryManaged<T> buffer3DT(refCoilMapN.get_dimensions(), gtPlus_mem_manager_);

                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(refCoilMapN, buffer3DT));

                hoNDArray<T> coilMapN(RO, E1, E2, CHA, workOrder3DT->coilMap_->begin()+usedN*RO*E1*E2*CHA);

                if ( workOrder3DT->csm_use_gpu_ )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIHGPU_FullResMap(buffer3DT, 
                            coilMapN, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(buffer3DT, 
                            coilMapN, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
                }
                GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder3DT->coilMap_, usedN));
            }
            else
            {
                hoNDArrayMemoryManaged<T> buffer3DT(workOrder3DT->ref_coil_map_.get_dimensions(), gtPlus_mem_manager_);

                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->ref_coil_map_, buffer3DT));

                if ( workOrder3DT->csm_use_gpu_ )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIHGPU_FullResMap(buffer3DT, 
                            *workOrder3DT->coilMap_, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(buffer3DT, 
                            *workOrder3DT->coilMap_, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
                }
            }

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *workOrder3DT->coilMap_, "coilMap_");
        }

        // partial fourier handling
        GADGET_CHECK_RETURN_FALSE(this->performPartialFourierHandling(workOrder3DT));

        workOrder3DT->complexIm_.create(RO, E1, E2, N);

        GADGET_CHECK_PERFORM(performTiming_, gt_timer1_.start("perform coil combination"));

        hoNDArrayMemoryManaged<T> buffer3DT(workOrder3DT->data_.get_dimensions(), gtPlus_mem_manager_);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->data_, buffer3DT);
        gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(buffer3DT, *workOrder3DT->coilMap_, workOrder3DT->complexIm_);

        GADGET_CHECK_PERFORM(performTiming_, gt_timer1_.stop());

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder3DT->complexIm_, "combined");
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker3DTNoAcceleration<T>::performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT) ... ");
        return false;
    }

    return true;
}

}}
