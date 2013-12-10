/** \file   gtPlusISMRMRDReconWorker2DTGRAPPA.h
    \brief  Implement the 2DT GRAPPA reconstruction
    \author Hui Xue
*/

#pragma once

#include "ismrmrd.h"
#include "GadgetronTimer.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorker2DT.h"
#include "gtPlusGRAPPA.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorker2DTGRAPPA : public gtPlusReconWorker2DT<T>
{
public:

    typedef gtPlusReconWorker2DT<T> BaseClass;

    gtPlusReconWorker2DTGRAPPA() : BaseClass() {}
    virtual ~gtPlusReconWorker2DTGRAPPA() {}

    virtual bool performCalibPrep(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, gtPlusReconWorkOrder2DT<T>* workOrder2DT);
    virtual bool performCalibImpl(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, gtPlusReconWorkOrder2DT<T>* workOrder2DT, size_t n, size_t usedS);

    virtual bool performUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& data);

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_mem_manager_;

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

    gtPlusGRAPPA<T> grappa_;
};

template <typename T> 
bool gtPlusReconWorker2DTGRAPPA<T>::
performCalibPrep(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, gtPlusReconWorkOrder2DT<T>* workOrder2DT)
{
    try
    {
        size_t RO = workOrder2DT->data_.get_size(0);
        size_t E1 = workOrder2DT->data_.get_size(1);
        size_t N = workOrder2DT->data_.get_size(3);
        size_t S = workOrder2DT->data_.get_size(4);

        size_t srcCHA = ref_src.get_size(2);

        size_t refRO = ref_dst.get_size(0);
        size_t refE1 = ref_dst.get_size(1);
        size_t refN = ref_dst.get_size(3);
        size_t dstCHA = ref_dst.get_size(2);

        std::vector<int> kE1, oE1;
        bool fitItself = true;
        GADGET_CHECK_RETURN_FALSE(grappa_.kerPattern(kE1, oE1, workOrder2DT->acceFactorE1_, workOrder2DT->grappa_kSize_E1_, fitItself));

        size_t kRO = workOrder2DT->grappa_kSize_RO_;
        size_t kNE1 = workOrder2DT->grappa_kSize_E1_;
        size_t oNE1 = oE1.size();

        workOrder2DT->kernel_->create(kRO, kNE1, srcCHA, dstCHA, oNE1, refN, S);
        workOrder2DT->kernelIm_->create(RO, E1, srcCHA, dstCHA, refN, S);
        workOrder2DT->unmixingCoeffIm_->create(RO, E1, srcCHA, refN, S);
        workOrder2DT->gfactor_.create(RO, E1, refN, S);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker2DTGRAPPA<T>::performCalibPrep(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DTGRAPPA<T>::
performCalibImpl(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, gtPlusReconWorkOrder2DT<T>* workOrder2DT, size_t n, size_t usedS)
{
    try
    {
        size_t RO = workOrder2DT->data_.get_size(0);
        size_t E1 = workOrder2DT->data_.get_size(1);
        size_t N = workOrder2DT->data_.get_size(3);
        size_t S = workOrder2DT->data_.get_size(4);

        size_t srcCHA = ref_src.get_size(2);

        size_t refRO = ref_dst.get_size(0);
        size_t refE1 = ref_dst.get_size(1);
        size_t refN = ref_dst.get_size(3);
        size_t dstCHA = ref_dst.get_size(2);

        std::vector<int> kE1, oE1;
        bool fitItself = true;
        GADGET_CHECK_RETURN_FALSE(grappa_.kerPattern(kE1, oE1, workOrder2DT->acceFactorE1_, workOrder2DT->grappa_kSize_E1_, fitItself));

        size_t kRO = workOrder2DT->grappa_kSize_RO_;
        size_t kNE1 = workOrder2DT->grappa_kSize_E1_;
        size_t oNE1 = oE1.size();

        ho3DArray<T> acsSrc(refRO, refE1, srcCHA, const_cast<T*>(ref_src.begin()+n*refRO*refE1*srcCHA+usedS*refRO*refE1*srcCHA*refN));
        ho3DArray<T> acsDst(refRO, refE1, dstCHA, const_cast<T*>(ref_dst.begin()+n*refRO*refE1*dstCHA+usedS*refRO*refE1*dstCHA*refN));

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, acsSrc, "acsSrc");
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, acsDst, "acsDst");

        ho5DArray<T> ker(kRO, kNE1, srcCHA, dstCHA, oNE1, workOrder2DT->kernel_->begin()+n*kRO*kNE1*srcCHA*dstCHA*oNE1+usedS*kRO*kNE1*srcCHA*dstCHA*oNE1*refN);
        grappa_.calib(acsSrc, acsDst, workOrder2DT->grappa_reg_lamda_, kRO, kE1, oE1, ker);

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, ker, "ker");

        hoNDArray<T> kIm(RO, E1, srcCHA, dstCHA, workOrder2DT->kernelIm_->begin()+n*RO*E1*srcCHA*dstCHA+usedS*RO*E1*srcCHA*dstCHA*refN);
        grappa_.imageDomainKernel(ker, kRO, kE1, oE1, RO, E1, kIm);

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kIm, "kIm");

        hoNDArray<T> coilMap(RO, E1, dstCHA, workOrder2DT->coilMap_->begin()+n*RO*E1*dstCHA+usedS*RO*E1*dstCHA*refN);
        hoNDArray<T> unmixC(RO, E1, srcCHA, workOrder2DT->unmixingCoeffIm_->begin()+n*RO*E1*srcCHA+usedS*RO*E1*srcCHA*refN);
        hoNDArray<T> gFactor(RO, E1, workOrder2DT->gfactor_.begin()+n*RO*E1+usedS*RO*E1*refN);

        this->unmixCoeff(kIm, coilMap, unmixC, gFactor);
        GADGET_CHECK_RETURN_FALSE(Gadgetron::scal(1.0/workOrder2DT->acceFactorE1_, gFactor));

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, unmixC, "unmixC");
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, gFactor, "gFactor");
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker2DTGRAPPA<T>::performCalibImpl(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DTGRAPPA<T>::
performUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& data_dst)
{
    try
    {
        int n;

        size_t RO = workOrder2DT->data_.get_size(0);
        size_t E1 = workOrder2DT->data_.get_size(1);
        size_t N = workOrder2DT->data_.get_size(3);
        size_t S = workOrder2DT->data_.get_size(4);

        size_t srcCHA = workOrder2DT->kernelIm_->get_size(2);
        size_t dstCHA = workOrder2DT->kernelIm_->get_size(3);

        size_t refN = workOrder2DT->kernelIm_->get_size(4);

        workOrder2DT->complexIm_.create(RO, E1, N, S);

        if ( workOrder2DT->downstream_coil_compression_ )
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(workOrder2DT->data_, buffer2DT_);
        }
        else
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(data_dst, buffer2DT_);
        }

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, buffer2DT_, "buffer2DT_");

        bool recon_kspace = false;

        if ( workOrder2DT->CalibMode_ == ISMRMRD_embedded )
        {
            if ( workOrder2DT->embedded_fullres_coilmap_ || workOrder2DT->embedded_ref_fillback_ )
            {
                recon_kspace = true;
            }
        }

        if ( workOrder2DT->CalibMode_ == ISMRMRD_separate )
        {
            if ( workOrder2DT->separate_fullres_coilmap_ )
            {
                recon_kspace = true;
            }
        }

        if ( workOrder2DT->recon_kspace_needed_ )
        {
            recon_kspace = true;
        }

        // if kspace is actually needed
        if ( recon_kspace )
        {
            workOrder2DT->fullkspace_ = data_dst;

            buffer2DT_unwrapping_.create(RO, E1, srcCHA, dstCHA);

            size_t usedS;
            for ( usedS=0; usedS<S; usedS++ )
            {
                if ( (refN<N) || (refN==1) )
                {
                    hoNDArray<T> kIm(RO, E1, srcCHA, dstCHA, workOrder2DT->kernelIm_->begin()+usedS*RO*E1*srcCHA*dstCHA*refN);
                    hoNDArray<T> aliasedIm(RO, E1, srcCHA, N, buffer2DT_.begin()+usedS*RO*E1*srcCHA*N);
                    hoNDArray<T> unwarppedIm(RO, E1, dstCHA, N, workOrder2DT->fullkspace_.begin()+usedS*RO*E1*dstCHA*N);

                    this->applyImageDomainKernelImage(aliasedIm, kIm, buffer2DT_unwrapping_, unwarppedIm);
                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, unwarppedIm, "unwarppedIm");
                }
                else
                {
                    #pragma omp parallel private(n)
                    {
                        hoNDArray<T> complexIm(RO, E1, dstCHA);

                        #pragma omp for
                        for ( n=0; n<(int)N; n++ )
                        {
                            hoNDArray<T> kIm(RO, E1, srcCHA, dstCHA, workOrder2DT->kernelIm_->begin()+n*RO*E1*srcCHA*dstCHA+usedS*RO*E1*srcCHA*dstCHA*refN);

                            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kIm, "kIm_n");

                            T* pIm2D = buffer2DT_.begin()+n*RO*E1*srcCHA+usedS*RO*E1*srcCHA*N;
                            hoNDArray<T> aliasedIm(RO, E1, srcCHA, pIm2D);

                            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, aliasedIm, "aliasedIm_n");

                            this->applyImageDomainKernelImage(aliasedIm, kIm, complexIm);
                            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, complexIm, "complexIm_n");

                            memcpy(workOrder2DT->fullkspace_.begin()+n*RO*E1*dstCHA+usedS*RO*E1*dstCHA*N, complexIm.begin(), sizeof(T)*RO*E1*dstCHA);
                        }
                    }
                }

                hoNDArray<T> unwarppedIm(RO, E1, dstCHA, N, workOrder2DT->fullkspace_.begin()+usedS*RO*E1*dstCHA*N);
                hoNDArray<T> combined(RO, E1, N, workOrder2DT->complexIm_.begin()+usedS*RO*E1*N);

                if ( refN == N )
                {
                    hoNDArray<T> coilMap(RO, E1, dstCHA, refN, workOrder2DT->coilMap_->begin()+usedS*RO*E1*dstCHA*refN);
                    gtPlusISMRMRDReconUtilComplex<T>().coilCombine(unwarppedIm, coilMap, combined);
                }
                else
                {
                    hoNDArray<T> coilMap(RO, E1, dstCHA, workOrder2DT->coilMap_->begin()+usedS*RO*E1*dstCHA*refN);
                    gtPlusISMRMRDReconUtilComplex<T>().coilCombine(unwarppedIm, coilMap, combined);
                }

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, combined, "combined");
            }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(workOrder2DT->fullkspace_);
        }
        else
        {
            size_t usedS;
            for ( usedS=0; usedS<S; usedS++ )
            {
                if ( (refN<N) || (refN==1) )
                {
                    hoNDArray<T> unmixCoeff(RO, E1, srcCHA, workOrder2DT->unmixingCoeffIm_->begin()+usedS*RO*E1*srcCHA*refN);
                    hoNDArray<T> aliasedIm(RO, E1, srcCHA, N, buffer2DT_.begin()+usedS*RO*E1*srcCHA*N);
                    hoNDArray<T> unwarppedIm(RO, E1, 1, N, workOrder2DT->complexIm_.begin()+usedS*RO*E1*N);

                    this->applyUnmixCoeffImage(aliasedIm, unmixCoeff, unwarppedIm);

                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, unwarppedIm, "unwarppedIm");
                }
                else
                {
                    // #pragma omp parallel for private(n)
                    for ( n=0; n<(int)N; n++ )
                    {
                        hoNDArray<T> unmixCoeff(RO, E1, srcCHA, workOrder2DT->unmixingCoeffIm_->begin()+n*RO*E1*srcCHA+usedS*RO*E1*srcCHA*refN);
                        hoNDArray<T> aliasedIm(RO, E1, srcCHA, buffer2DT_.begin()+n*RO*E1*srcCHA+usedS*RO*E1*srcCHA*N);
                        hoNDArray<T> unwarppedIm(RO, E1, 1, workOrder2DT->complexIm_.begin()+n*RO*E1+usedS*RO*E1*N);

                        this->applyUnmixCoeffImage(aliasedIm, unmixCoeff, unwarppedIm);

                        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, unwarppedIm, "unwarppedIm");
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker2DTGRAPPA<T>::performUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& data) ... ");
        return false;
    }

    return true;
}

}}
