/** \file   gtPlusISMRMRDReconWorker3DTGRAPPA.h
    \brief  Implement the 3DT GRAPPA reconstruction
    \author Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"
#include "GadgetronTimer.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorker3DT.h"
#include "gtPlusGRAPPA.h"
#include "mri_core_grappa.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorker3DTGRAPPA : public gtPlusReconWorker3DT<T>
{
public:

    typedef gtPlusReconWorker3DT<T> BaseClass;
    typedef typename BaseClass::value_type value_type;

    typedef gtPlusReconWorkOrder3DT<T> WorkOrderType;

    gtPlusReconWorker3DTGRAPPA() : BaseClass() {}
    virtual ~gtPlusReconWorker3DTGRAPPA() {}

    virtual bool performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT);

    virtual bool performCalibPrep(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, WorkOrderType* workOrder3DT);
    virtual bool performCalibImpl(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, WorkOrderType* workOrder3DT, size_t usedN);

    virtual bool performUnwrapping(gtPlusReconWorkOrder3DT<T>* workOrder3DT, const hoNDArray<T>& data);

    virtual bool computeKSpace(gtPlusReconWorkOrder3DT<T>* workOrder3DT);

    virtual bool splitJob(gtPlusReconWorkOrder3DT<T>* workOrder3DT, size_t& jobN);

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

    gtPlusGRAPPA<T> grappa_;
};

template <typename T> 
bool gtPlusReconWorker3DTGRAPPA<T>::computeKSpace(gtPlusReconWorkOrder3DT<T>* workOrder3DT)
{
    bool recon_kspace = false;

    if ( workOrder3DT->CalibMode_ == ISMRMRD_embedded )
    {
        if ( workOrder3DT->embedded_fullres_coilmap_ || workOrder3DT->embedded_ref_fillback_ )
        {
            recon_kspace = true;
        }
    }

    if ( workOrder3DT->CalibMode_ == ISMRMRD_separate )
    {
        if ( workOrder3DT->separate_fullres_coilmap_ )
        {
            recon_kspace = true;
        }
    }

    if ( workOrder3DT->recon_kspace_needed_ )
    {
        recon_kspace = true;
    }

    return recon_kspace;
}

template <typename T> 
bool gtPlusReconWorker3DTGRAPPA<T>::
splitJob(gtPlusReconWorkOrder3DT<T>* workOrder3DT, size_t& jobN)
{
    return false;

    size_t RO = workOrder3DT->data_.get_size(0);
    size_t E1 = workOrder3DT->data_.get_size(1);
    size_t E2 = workOrder3DT->data_.get_size(2);

    size_t srcCHA = workOrder3DT->kernel_->get_size(3);
    size_t dstCHA = workOrder3DT->kernel_->get_size(4);

    jobN = workOrder3DT->job_num_of_N_;
    size_t jobMegaBytes = workOrder3DT->job_max_Megabytes_;

    bool splitJobs = (jobN>0 && RO>jobN);
    if ( !splitJobs )
    {
        if ( jobMegaBytes>0 )
        {
            size_t jobN = jobMegaBytes/(E1*E2*srcCHA*dstCHA*sizeof(T)/1024/1024);
            if ( jobN < RO ) splitJobs = true;
            GDEBUG_STREAM("grappa - 3DT - size of largest job : " << jobN);
        }
    }

    bool reconKSpace = this->computeKSpace(workOrder3DT);
    if ( !reconKSpace )
    {
        splitJobs = false;
    }

    return splitJobs;
}

template <typename T> 
bool gtPlusReconWorker3DTGRAPPA<T>::
performCalibPrep(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, WorkOrderType* workOrder3DT)
{
    grappa_.performTiming_ = performTiming_;

    size_t RO = workOrder3DT->data_.get_size(0);
    size_t E1 = workOrder3DT->data_.get_size(1);
    size_t E2 = workOrder3DT->data_.get_size(2);
    size_t N = workOrder3DT->data_.get_size(4);
    size_t srcCHA = workOrder3DT->data_.get_size(3);

    size_t refRO = ref_dst.get_size(0);
    size_t refE1 = ref_dst.get_size(1);
    size_t refE2 = ref_dst.get_size(2);
    size_t refN = ref_dst.get_size(4);
    size_t dstCHA = ref_dst.get_size(3);

    bool reconKSpace = this->computeKSpace(workOrder3DT);

    std::vector<int> kE1, oE1;
    bool fitItself = true;
    // GADGET_CHECK_RETURN_FALSE(grappa_.kerPattern(kE1, oE1, (int)workOrder3DT->acceFactorE1_, workOrder3DT->grappa_kSize_E1_, fitItself));

    std::vector<int> kE2, oE2;
    // GADGET_CHECK_RETURN_FALSE(grappa_.kerPattern(kE2, oE2, (int)workOrder3DT->acceFactorE2_, workOrder3DT->grappa_kSize_E2_, fitItself));

    size_t convKRO, convKE1, convKE2;

    size_t kRO = workOrder3DT->grappa_kSize_RO_;
    size_t kNE1 = workOrder3DT->grappa_kSize_E1_;
    size_t kNE2 = workOrder3DT->grappa_kSize_E2_;

    grappa3d_kerPattern(kE1, oE1, kE2, oE2, convKRO, convKE1, convKE2,
        (int)workOrder3DT->acceFactorE1_, workOrder3DT->acceFactorE2_,
        kRO, kNE1, kNE2, fitItself);

    workOrder3DT->kernel_->create(convKRO, convKE1, convKE2, srcCHA, dstCHA, refN);
    Gadgetron::clear(workOrder3DT->kernel_.get());

    /*if ( performTiming_ ) { gt_timer3_.start("allocate image domain kernel ... "); }
    workOrder3DT->kernelIm_->create(RO, E1, E2, srcCHA, dstCHA, refN);
    if ( performTiming_ ) { gt_timer3_.stop(); }*/

    if ( !reconKSpace )
    {
        if ( performTiming_ ) { gt_timer3_.start("allocate unmixing coefficient ... "); }
        workOrder3DT->unmixingCoeffIm_->create(RO, E1, E2, srcCHA, refN);
        Gadgetron::clear(workOrder3DT->unmixingCoeffIm_.get());
        if ( performTiming_ ) { gt_timer3_.stop(); }

        workOrder3DT->gfactor_.create(RO, E1, E2, refN);
        Gadgetron::clear(&(workOrder3DT->gfactor_));
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTGRAPPA<T>::
performCalibImpl(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, WorkOrderType* workOrder3DT, size_t usedN)
{
    size_t RO = workOrder3DT->data_.get_size(0);
    size_t E1 = workOrder3DT->data_.get_size(1);
    size_t E2 = workOrder3DT->data_.get_size(2);
    size_t N = workOrder3DT->data_.get_size(4);
    size_t srcCHA = workOrder3DT->data_.get_size(3);

    size_t refRO = ref_dst.get_size(0);
    size_t refE1 = ref_dst.get_size(1);
    size_t refE2 = ref_dst.get_size(2);
    size_t refN = ref_dst.get_size(4);
    size_t dstCHA = ref_dst.get_size(3);

    bool reconKSpace = this->computeKSpace(workOrder3DT);

    ho4DArray<T> acsSrc(refRO, refE1, refE2, srcCHA, const_cast<T*>(ref_src.begin()+usedN*refRO*refE1*refE2*srcCHA));
    ho4DArray<T> acsDst(refRO, refE1, refE2, dstCHA, const_cast<T*>(ref_dst.begin()+usedN*refRO*refE1*refE2*dstCHA));

    std::ostringstream ostr;
    ostr << "_n_" << usedN;
    std::string suffix = ostr.str();

    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(acsSrc, debugFolder_+"acsSrc"+suffix); }
    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(acsDst, debugFolder_+"acsDst"+suffix); }

    size_t kRO = workOrder3DT->grappa_kSize_RO_;
    size_t kNE1 = workOrder3DT->grappa_kSize_E1_;
    size_t kNE2 = workOrder3DT->grappa_kSize_E2_;

    size_t convKRO = workOrder3DT->kernel_->get_size(0);
    size_t convKNE1 = workOrder3DT->kernel_->get_size(1);
    size_t convKNE2 = workOrder3DT->kernel_->get_size(2);

    hoNDArray<T> ker(convKRO, convKNE1, convKNE2, srcCHA, dstCHA, workOrder3DT->kernel_->begin() + usedN*convKRO*convKNE1*convKNE2*srcCHA*dstCHA);

    if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration ... "); }
    Gadgetron::grappa3d_calib_convolution_kernel(acsSrc, acsDst, workOrder3DT->acceFactorE1_, workOrder3DT->acceFactorE2_, workOrder3DT->grappa_reg_lamda_, workOrder3DT->grappa_calib_over_determine_ratio_, kRO, kNE1, kNE2, ker);
    if ( performTiming_ ) { gt_timer3_.stop(); }

    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ker, debugFolder_+"ker"+suffix); }

    if ( !reconKSpace )
    {
        hoNDArray<T> coilMap(RO, E1, E2, dstCHA, workOrder3DT->coilMap_->begin()+usedN*RO*E1*E2*dstCHA);
        hoNDArray<T> unmixC(RO, E1, E2, srcCHA);
        hoNDArray<T> gFactor(RO, E1, E2, 1, workOrder3DT->gfactor_.begin()+usedN*RO*E1*E2);

        //this->unmixCoeff(kIm, coilMap, unmixC, gFactor);
        //GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (value_type)(1.0/workOrder3DT->acceFactorE1_/workOrder3DT->acceFactorE2_), gFactor));

        hoNDArray<value_type> gFactorReal(RO, E1, E2);

        if (performTiming_) { gt_timer3_.start("grappa 3D, compute unmixing coefficient ... "); }
        Gadgetron::grappa3d_unmixing_coeff(ker, coilMap, workOrder3DT->acceFactorE1_, workOrder3DT->acceFactorE2_, unmixC, gFactorReal);
        if (performTiming_) { gt_timer3_.stop(); }

        Gadgetron::real_to_complex(gFactorReal, gFactor);

        memcpy(workOrder3DT->unmixingCoeffIm_->begin()+usedN*RO*E1*E2*srcCHA, unmixC.begin(), unmixC.get_number_of_bytes());

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(unmixC, debugFolder_+"unmixC"+suffix); }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(gFactor, debugFolder_+"gFactor"+suffix); }
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTGRAPPA<T>::
performUnwrapping(gtPlusReconWorkOrder3DT<T>* workOrder3DT, const hoNDArray<T>& data_dst)
{
    try
    {
        int e1, e2, n;

        size_t RO = workOrder3DT->data_.get_size(0);
        size_t E1 = workOrder3DT->data_.get_size(1);
        size_t E2 = workOrder3DT->data_.get_size(2);
        size_t N = workOrder3DT->data_.get_size(4);

        size_t kRO = workOrder3DT->kernel_->get_size(0);
        size_t kNE1 = workOrder3DT->kernel_->get_size(1);
        size_t kNE2 = workOrder3DT->kernel_->get_size(2);

        size_t srcCHA = workOrder3DT->kernel_->get_size(3);
        size_t dstCHA = workOrder3DT->kernel_->get_size(4);

        size_t refN = workOrder3DT->kernel_->get_size(5);

        //size_t srcCHA = workOrder3DT->kernelIm_->get_size(3);
        //size_t dstCHA = workOrder3DT->kernelIm_->get_size(4);

        //size_t refN = workOrder3DT->kernelIm_->get_size(5);

        workOrder3DT->complexIm_.create(RO, E1, E2, 1, N);

        hoNDArray<T> aliasedIm;

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D compute aliased image ... "); }
        if ( workOrder3DT->downstream_coil_compression_ )
        {
            aliasedIm.create(workOrder3DT->data_.get_dimensions());
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->data_, aliasedIm);
        }
        else
        {
            aliasedIm.create(data_dst.get_dimensions());
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(data_dst, aliasedIm);
        }

        // find the effective acce factor
        size_t num = 0;
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
            if (this->verbose_) GDEBUG_STREAM("gtPlusReconWorker3DTGRAPPA, length for RO : " << lenRO << " - " << lenRO/RO);

            double effectiveAcceFactor = (double)(N*E1*E2) / (num);
            if (this->verbose_) GDEBUG_STREAM("gtPlusReconWorker3DTGRAPPA, effectiveAcceFactor : " << effectiveAcceFactor);

            double ROScalingFactor = (double)RO / (double)lenRO;

            // since the grappa in gadgetron is doing signal preserving scaling, to perserve noise level, we need this compensation factor
            double grappaKernelCompensationFactor = 1.0 / (workOrder3DT->acceFactorE1_*workOrder3DT->acceFactorE2_);

            typename realType<T>::Type fftCompensationRatio = (typename realType<T>::Type)(std::sqrt(ROScalingFactor*effectiveAcceFactor) * grappaKernelCompensationFactor);

            if (this->verbose_) GDEBUG_STREAM("gtPlusReconWorker3DTGRAPPA, fftCompensationRatio : " << fftCompensationRatio);

            Gadgetron::scal(fftCompensationRatio, aliasedIm);

            // if the image data is scaled and ref lines are going to be filled back to the data, 
            // the reference lines should be scaled too
            if (workOrder3DT->CalibMode_ == ISMRMRD_embedded)
            {
                if (workOrder3DT->embedded_ref_fillback_)
                {
                    if (this->verbose_) GDEBUG_STREAM("gtPlusReconWorker3DTGRAPPA, ref fill back, scaling ref : " << fftCompensationRatio);
                    Gadgetron::scal(fftCompensationRatio, workOrder3DT->ref_);
                }
            }
        }

        //typename realType<T>::Type fftCompensationRatio = (typename realType<T>::Type)(1.0/std::sqrt( (double)workOrder3DT->acceFactorE1_ * (double)workOrder3DT->acceFactorE2_ ));
        //Gadgetron::scal( fftCompensationRatio, aliasedIm);

        //// if the image data is scaled and ref lines are going to be filled back to the data, 
        //// the reference lines should be scaled too
        //if ( workOrder3DT->CalibMode_ == ISMRMRD_embedded )
        //{
        //    if ( workOrder3DT->embedded_ref_fillback_ )
        //    {
        //        Gadgetron::scal( fftCompensationRatio, workOrder3DT->ref_);
        //    }
        //}

        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(aliasedIm, debugFolder_+"aliasedIm"); }

        bool recon_kspace = this->computeKSpace(workOrder3DT);

        // if kspace is actually needed
        if ( recon_kspace )
        {
            workOrder3DT->fullkspace_ = data_dst;

            if ( (refN<N) || (refN==1) )
            {
                hoNDArray<T> convKer(kRO, kNE1, kNE2, srcCHA, dstCHA, workOrder3DT->kernel_->begin());

                if ( performTiming_ ) { gt_timer3_.start("grappa 3D apply image domain kernel for every channel ... "); }
                Gadgetron::grappa3d_image_domain_unwrapping_aliasedImage(convKer, aliasedIm, workOrder3DT->acceFactorE1_, workOrder3DT->acceFactorE2_, workOrder3DT->fullkspace_);
                if ( performTiming_ ) { gt_timer3_.stop(); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->fullkspace_, debugFolder_+"unwarppedIm"); }
            }
            else
            {
                hoNDArray<T> complexIm(RO, E1, E2, dstCHA);
                for ( n=0; n<(int)N; n++ )
                {
                    hoNDArray<T> convKer(kRO, kNE1, kNE2, srcCHA, dstCHA, workOrder3DT->kernel_->begin() + n*kRO*kNE1*kNE2*srcCHA*dstCHA);

                    if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(convKer, debugFolder_ + "convKer_n"); }

                    hoNDArray<T> aliasedImN(RO, E1, E2, srcCHA, aliasedIm.begin()+n*RO*E1*E2*srcCHA);

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(aliasedImN, debugFolder_+"aliasedIm_n"); }

                    Gadgetron::grappa3d_image_domain_unwrapping_aliasedImage(convKer, aliasedImN, workOrder3DT->acceFactorE1_, workOrder3DT->acceFactorE2_, complexIm);
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(complexIm, debugFolder_+"complexIm_n"); }

                    memcpy(workOrder3DT->fullkspace_.begin()+n*RO*E1*E2*dstCHA, complexIm.begin(), sizeof(T)*RO*E1*E2*dstCHA);
                }
            }

            if ( (workOrder3DT->coilMap_->get_size(0)==RO) 
                && (workOrder3DT->coilMap_->get_size(1)==E1) 
                && (workOrder3DT->coilMap_->get_size(2)==E2) 
                && (workOrder3DT->coilMap_->get_size(3)==dstCHA) )
            {
                if ( performTiming_ ) { gt_timer3_.start("grappa 3D coil combination ... "); }
                // gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(workOrder3DT->fullkspace_, *workOrder3DT->coilMap_, workOrder3DT->complexIm_);
                Gadgetron::coil_combine(workOrder3DT->fullkspace_, *workOrder3DT->coilMap_, 3, workOrder3DT->complexIm_);
                if ( performTiming_ ) { gt_timer3_.stop(); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"combined"); }
            }

            if ( performTiming_ ) { gt_timer3_.start("grappa 3D go back to kspace ... "); }
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(workOrder3DT->fullkspace_);
            if ( performTiming_ ) { gt_timer3_.stop(); }
        }
        else
        {
            if ( (refN<N) || (refN==1) )
            {
                hoNDArray<T> unmixCoeff(RO, E1, E2, srcCHA, workOrder3DT->unmixingCoeffIm_->begin());

                if ( performTiming_ ) { gt_timer3_.start("grappa 3D apply unmixing coeff ... "); }
                Gadgetron::apply_unmix_coeff_aliased_image_3D(aliasedIm, unmixCoeff, workOrder3DT->complexIm_);
                if ( performTiming_ ) { gt_timer3_.stop(); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"unwarppedIm"); }
            }
            else
            {
                for ( n=0; n<(int)N; n++ )
                {
                    hoNDArray<T> unmixCoeff(RO, E1, E2, srcCHA, workOrder3DT->unmixingCoeffIm_->begin()+n*RO*E1*E2*srcCHA);
                    hoNDArray<T> aliasedImN(RO, E1, E2, srcCHA, aliasedIm.begin()+n*RO*E1*E2*srcCHA);
                    hoNDArray<T> unwarppedIm(RO, E1, E2, 1, workOrder3DT->complexIm_.begin()+n*RO*E1*E2);

                    Gadgetron::apply_unmix_coeff_aliased_image_3D(aliasedImN, unmixCoeff, unwarppedIm);

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(unwarppedIm, debugFolder_+"unwarppedIm"); }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTGRAPPA<T>::performUnwrapping(gtPlusReconWorkOrder3DT<T>* workOrder3DT, const hoNDArray<T>& data) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTGRAPPA<T>::performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(workOrder3DT!=NULL);

        // call the BaseClass
        GADGET_CHECK_RETURN_FALSE(BaseClass::performRecon(workOrder3DT));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTGRAPPA<T>::performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT) ... ");
        return false;
    }

    return true;
}

}}
