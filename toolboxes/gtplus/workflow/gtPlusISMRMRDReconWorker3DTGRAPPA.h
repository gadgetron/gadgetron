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
    using BaseClass::gtPlus_mem_manager_;

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
    GADGET_CHECK_RETURN_FALSE(grappa_.kerPattern(kE1, oE1, (int)workOrder3DT->acceFactorE1_, workOrder3DT->grappa_kSize_E1_, fitItself));

    std::vector<int> kE2, oE2;
    GADGET_CHECK_RETURN_FALSE(grappa_.kerPattern(kE2, oE2, (int)workOrder3DT->acceFactorE2_, workOrder3DT->grappa_kSize_E2_, fitItself));

    size_t kRO = workOrder3DT->grappa_kSize_RO_;
    size_t kNE1 = workOrder3DT->grappa_kSize_E1_;
    size_t oNE1 = oE1.size();

    size_t kNE2 = workOrder3DT->grappa_kSize_E2_;
    size_t oNE2 = oE1.size();

    workOrder3DT->kernel_->create(kRO, kNE1, kNE2, srcCHA, dstCHA, oNE1, oNE2, refN);
    Gadgetron::clear(workOrder3DT->kernel_.get());

    size_t jobN;
    bool splitJobs = this->splitJob(workOrder3DT, jobN);

    if ( !splitJobs )
    {
        if ( performTiming_ ) { gt_timer3_.start("allocate image domain kernel ... "); }
        if ( gtPlus_mem_manager_ )
        {
            if ( workOrder3DT->kernelIm_->get_number_of_elements() != (size_t)RO*E1*E2*srcCHA*dstCHA*refN )
            {
                workOrder3DT->kernelIm_->create(RO, E1, E2, srcCHA, dstCHA, refN, (T*)(gtPlus_mem_manager_->allocate(sizeof(T)*(size_t)RO*E1*E2*srcCHA*dstCHA*refN)));
            }
        }
        else
        {
            workOrder3DT->kernelIm_->create(RO, E1, E2, srcCHA, dstCHA, refN);
        }
        if ( performTiming_ ) { gt_timer3_.stop(); }
    }
    else
    {
        int maxKE1 = std::abs(kE1[0]);
        if ( std::abs(kE1[kNE1-1]) > maxKE1 )
        {
            maxKE1 = std::abs(kE1[kNE1-1]);
        }
        int convKE1 = 2*maxKE1+1;

        int maxKE2 = std::abs(kE2[0]);
        if ( std::abs(kE2[kNE2-1]) > maxKE2 )
        {
            maxKE2 = std::abs(kE2[kNE2-1]);
        }
        int convKE2 = 2*maxKE2+1;

        if ( performTiming_ ) { gt_timer3_.start("allocate image domain kernel only along RO ... "); }
        if ( gtPlus_mem_manager_ )
        {
            if ( workOrder3DT->kernelIm_->get_number_of_elements() != (size_t)RO*convKE1*convKE2*srcCHA*dstCHA*refN )
            {
                workOrder3DT->kernelIm_->create(convKE1, convKE2, RO, srcCHA, dstCHA, refN, (T*)(gtPlus_mem_manager_->allocate(sizeof(T)*(size_t)RO*convKE1*convKE2*srcCHA*dstCHA*refN)));
            }
        }
        else
        {
            workOrder3DT->kernelIm_->create(convKE1, convKE2, RO, srcCHA, dstCHA, refN);
            // pre-set to zero is needed here
            memset(workOrder3DT->kernelIm_->begin(), 0, workOrder3DT->kernelIm_->get_number_of_bytes());
        }
        if ( performTiming_ ) { gt_timer3_.stop(); }
    }

    if ( !reconKSpace )
    {
        if ( performTiming_ ) { gt_timer3_.start("allocate unmixing coefficient ... "); }
        if ( gtPlus_mem_manager_ )
        {
            if ( workOrder3DT->unmixingCoeffIm_->get_number_of_elements() != (size_t)RO*E1*E2*srcCHA*refN )
            {
                workOrder3DT->unmixingCoeffIm_->create(RO, E1, E2, srcCHA, refN, (T*)(gtPlus_mem_manager_->allocate(sizeof(T)*(size_t)RO*E1*E2*srcCHA*refN)));
            }
        }
        else
        {
            workOrder3DT->unmixingCoeffIm_->create(RO, E1, E2, srcCHA, refN);
        }
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

    std::vector<int> kE1, oE1;
    bool fitItself = true;
    GADGET_CHECK_RETURN_FALSE(grappa_.kerPattern(kE1, oE1, (size_t)workOrder3DT->acceFactorE1_, workOrder3DT->grappa_kSize_E1_, fitItself));

    std::vector<int> kE2, oE2;
    GADGET_CHECK_RETURN_FALSE(grappa_.kerPattern(kE2, oE2, (size_t)workOrder3DT->acceFactorE2_, workOrder3DT->grappa_kSize_E2_, fitItself));

    size_t kRO = workOrder3DT->grappa_kSize_RO_;
    size_t kNE1 = workOrder3DT->grappa_kSize_E1_;
    size_t oNE1 = oE1.size();

    size_t kNE2 = workOrder3DT->grappa_kSize_E2_;
    size_t oNE2 = oE1.size();

    ho4DArray<T> acsSrc(refRO, refE1, refE2, srcCHA, const_cast<T*>(ref_src.begin()+usedN*refRO*refE1*refE2*srcCHA));
    ho4DArray<T> acsDst(refRO, refE1, refE2, dstCHA, const_cast<T*>(ref_dst.begin()+usedN*refRO*refE1*refE2*dstCHA));

    std::ostringstream ostr;
    ostr << "_n_" << usedN;
    std::string suffix = ostr.str();

    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, acsSrc, "acsSrc"+suffix);
    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, acsDst, "acsDst"+suffix);

    grappa_.calib_use_gpu_  = workOrder3DT->grappa_use_gpu_;

    ho7DArray<T> ker(kRO, kNE1, kNE2, srcCHA, dstCHA, oNE1, oNE2, workOrder3DT->kernel_->begin()+usedN*kRO*kNE1*kNE2*srcCHA*dstCHA*oNE1*oNE2);
    if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration ... "); }
    grappa_.calib3D(acsSrc, acsDst, workOrder3DT->grappa_reg_lamda_, workOrder3DT->grappa_calib_over_determine_ratio_, kRO, kE1, kE2, oE1, oE2, ker);
    if ( performTiming_ ) { gt_timer3_.stop(); }

    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, ker, "ker"+suffix);

    size_t jobN;
    bool splitJobs = this->splitJob(workOrder3DT, jobN);

    if ( !splitJobs )
    {
        hoNDArray<T> kIm(RO, E1, E2, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin()+usedN*RO*E1*E2*srcCHA*dstCHA);
        if ( performTiming_ ) { gt_timer3_.start("grappa 3D image domain kernel ... "); }
        grappa_.imageDomainKernel3D(ker, kRO, kE1, kE2, oE1, oE2, RO, E1, E2, kIm);
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !reconKSpace )
        {
            hoNDArray<T> coilMap(RO, E1, E2, dstCHA, workOrder3DT->coilMap_->begin()+usedN*RO*E1*E2*dstCHA);
            hoNDArray<T> unmixC(RO, E1, E2, srcCHA);
            hoNDArray<T> gFactor(RO, E1, E2, workOrder3DT->gfactor_.begin()+usedN*RO*E1*E2);

            this->unmixCoeff(kIm, coilMap, unmixC, gFactor);
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (value_type)(1.0/workOrder3DT->acceFactorE1_/workOrder3DT->acceFactorE2_), gFactor));

            memcpy(workOrder3DT->unmixingCoeffIm_->begin()+usedN*RO*E1*E2*srcCHA, unmixC.begin(), unmixC.get_number_of_bytes());

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, unmixC, "unmixC"+suffix);
            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, gFactor, "gFactor"+suffix);
        }
    }
    else
    {
        int maxKE1 = std::abs(kE1[0]);
        if ( std::abs(kE1[kNE1-1]) > maxKE1 )
        {
            maxKE1 = std::abs(kE1[kNE1-1]);
        }
        int convKE1 = 2*maxKE1+1;

        int maxKE2 = std::abs(kE2[0]);
        if ( std::abs(kE2[kNE2-1]) > maxKE2 )
        {
            maxKE2 = std::abs(kE2[kNE2-1]);
        }
        int convKE2 = 2*maxKE2+1;

        hoNDArray<T> kIm(convKE1, convKE2, RO, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin()+usedN*convKE1*convKE2*RO*srcCHA*dstCHA);

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D image domain kernel only along RO ... "); }
        GADGET_CHECK_RETURN_FALSE(grappa_.imageDomainKernelRO3D(ker, kRO, kE1, kE2, oE1, oE2, RO, kIm));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() )
        {
            hoNDArray<T> kImROACha(convKE1, convKE2, RO, srcCHA, kIm.begin());
            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kImROACha, "kImROACha"+suffix);
        }
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTGRAPPA<T>::
performUnwrapping(gtPlusReconWorkOrder3DT<T>* workOrder3DT, const hoNDArray<T>& data_dst)
{
    try
    {
        int n;

        size_t RO = workOrder3DT->data_.get_size(0);
        size_t E1 = workOrder3DT->data_.get_size(1);
        size_t E2 = workOrder3DT->data_.get_size(2);
        size_t N = workOrder3DT->data_.get_size(4);

        size_t srcCHA = workOrder3DT->kernelIm_->get_size(3);
        size_t dstCHA = workOrder3DT->kernelIm_->get_size(4);

        size_t refN = workOrder3DT->kernelIm_->get_size(5);

        workOrder3DT->complexIm_.create(RO, E1, E2, 1, N);

        hoNDArrayMemoryManaged<T> aliasedIm(gtPlus_mem_manager_);

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

        typename realType<T>::Type fftCompensationRatio = (typename realType<T>::Type)(1.0/std::sqrt( (double)workOrder3DT->acceFactorE1_ * (double)workOrder3DT->acceFactorE2_ ));
        Gadgetron::scal( fftCompensationRatio, aliasedIm);

        // if the image data is scaled and ref lines are going to be filled back to the data, 
        // the reference lines should be scaled too
        if ( workOrder3DT->CalibMode_ == ISMRMRD_embedded )
        {
            if ( workOrder3DT->embedded_ref_fillback_ )
            {
                Gadgetron::scal( fftCompensationRatio, workOrder3DT->ref_);
            }
        }

        if ( performTiming_ ) { gt_timer3_.stop(); }

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, aliasedIm, "aliasedIm");

        bool recon_kspace = this->computeKSpace(workOrder3DT);

        // if kspace is actually needed
        if ( recon_kspace )
        {
            workOrder3DT->fullkspace_ = data_dst;

            size_t jobN;
            bool splitJobs = this->splitJob(workOrder3DT, jobN);

            if ( splitJobs )
            {
                size_t kE1 = workOrder3DT->kernelIm_->get_size(0);
                size_t kE2 = workOrder3DT->kernelIm_->get_size(1);
                size_t kRO = workOrder3DT->kernelIm_->get_size(2);

                if ( (refN<N) || (refN==1) )
                {
                    hoNDArray<T> kImPermuted(kE1, kE2, RO, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin());

                    hoNDArray<T> kImPermutedJob(kE1, kE2, jobN, srcCHA, dstCHA);

                    if ( performTiming_ ) { gt_timer3_.start("grappa 3D allocate buffer for kImPermutedZeroFilledJob ... "); }
                    hoNDArrayMemoryManaged<T> kImPermutedZeroFilledJob(E1, E2, jobN, srcCHA, dstCHA, gtPlus_mem_manager_);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    // aliased images
                    hoNDArray<T> aliasedImPermutedJob(E1, E2, jobN, srcCHA);

                    if ( performTiming_ ) { gt_timer3_.start("grappa 3D allocate buffer for aliasedIm permuted ... "); }
                    hoNDArrayMemoryManaged<T> aliasedImPermuted(E1, E2, RO, srcCHA, N, gtPlus_mem_manager_);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( performTiming_ ) { gt_timer3_.start("permuteROTo3rdDimensionFor3DRecon for aliased images ... "); }
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo3rdDimensionFor3DRecon(aliasedIm, aliasedImPermuted));
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    // unwrapped images
                    hoNDArray<T> unwrappedImPermutedJob(E1, E2, jobN, srcCHA, N);

                    if ( performTiming_ ) { gt_timer3_.start("grappa 3D allocate buffer for unwrapped images permuted ... "); }
                    hoNDArrayMemoryManaged<T> unwrappedImPermuted(E1, E2, RO, dstCHA, N, gtPlus_mem_manager_);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    // buffer
                    if ( performTiming_ ) { gt_timer3_.start("grappa 3D allocate buffer for unwrapping ... "); }
                    hoNDArrayMemoryManaged<T> buffer3DT_unwrapping(E1, E2, jobN, srcCHA, gtPlus_mem_manager_);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    size_t ro=0;
                    while ( ro<RO )
                    {
                        size_t start = ro;
                        size_t end = ro+jobN-1;
                        if ( end >= RO )
                        {
                            end = RO-1;
                            start = end-jobN+1;
                        }

                        GDEBUG_STREAM("grappa 3D - processing " << start << " to " << end << " ... ");

                        if ( (refN<N) || (refN==1) )
                        {
                            hoNDArray<T> kImPermuted(kE1, kE2, RO, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin());

                            if ( performTiming_ ) { gt_timer3_.start("cropOver3rdDimension hybrid domain kernel ... "); }
                            GADGET_CHECK_RETURN_FALSE(cropOver3rdDimension(kImPermuted, kImPermutedJob, start, end));
                            if ( performTiming_ ) { gt_timer3_.stop(); }

                            if ( performTiming_ ) { gt_timer3_.start("imageDomainKernelE1E2RO ... "); }
                            GADGET_CHECK_RETURN_FALSE(grappa_.imageDomainKernelE1E2RO(kImPermutedJob, (int)E1, (int)E2, kImPermutedZeroFilledJob));
                            if ( performTiming_ ) { gt_timer3_.stop(); }

                            if ( performTiming_ ) { gt_timer3_.start("cropOver3rdDimension aliased images ... "); }
                            GADGET_CHECK_RETURN_FALSE(cropOver3rdDimension(aliasedImPermuted, aliasedImPermutedJob, start, end));
                            if ( performTiming_ ) { gt_timer3_.stop(); }

                            if ( performTiming_ ) { gt_timer3_.start("grappa 3D apply image domain kernel for every channel and every job ... "); }
                            this->applyImageDomainKernelImage(aliasedImPermutedJob, kImPermutedZeroFilledJob, buffer3DT_unwrapping, unwrappedImPermutedJob);
                            if ( performTiming_ ) { gt_timer3_.stop(); }

                            if ( performTiming_ ) { gt_timer3_.start("setSubArrayOver3rdDimension unwrapped images ... "); }
                            GADGET_CHECK_RETURN_FALSE(setSubArrayOver3rdDimension(unwrappedImPermutedJob, unwrappedImPermuted, start, end));
                            if ( performTiming_ ) { gt_timer3_.stop(); }
                        }
                        else
                        {
                            for ( n=0; n<(int)N; n++ )
                            {
                                hoNDArray<T> kImPermuted(kE1, kE2, RO, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin()+n*kE1*kE2*RO*srcCHA*dstCHA);

                                if ( performTiming_ ) { gt_timer3_.start("cropOver3rdDimension hybrid domain kernel ... "); }
                                GADGET_CHECK_RETURN_FALSE(cropOver3rdDimension(kImPermuted, kImPermutedJob, start, end));
                                if ( performTiming_ ) { gt_timer3_.stop(); }

                                if ( performTiming_ ) { gt_timer3_.start("imageDomainKernelE1E2RO ... "); }
                                GADGET_CHECK_RETURN_FALSE(grappa_.imageDomainKernelE1E2RO(kImPermutedJob, (int)E1, (int)E2, kImPermutedZeroFilledJob));
                                if ( performTiming_ ) { gt_timer3_.stop(); }

                                hoNDArray<T> aliasedImPermutedN(E1, E2, RO, srcCHA, aliasedImPermuted.begin()+n*E1*E2*RO*srcCHA);

                                if ( performTiming_ ) { gt_timer3_.start("cropOver3rdDimension aliased images ... "); }
                                GADGET_CHECK_RETURN_FALSE(cropOver3rdDimension(aliasedImPermutedN, aliasedImPermutedJob, start, end));
                                if ( performTiming_ ) { gt_timer3_.stop(); }

                                if ( performTiming_ ) { gt_timer3_.start("grappa 3D apply image domain kernel for every channel and every job ... "); }
                                this->applyImageDomainKernelImage(aliasedImPermutedJob, kImPermutedZeroFilledJob, buffer3DT_unwrapping, unwrappedImPermutedJob);
                                if ( performTiming_ ) { gt_timer3_.stop(); }

                                if ( performTiming_ ) { gt_timer3_.start("setSubArrayOver3rdDimension unwrapped images ... "); }
                                GADGET_CHECK_RETURN_FALSE(setSubArrayOver3rdDimension(unwrappedImPermutedJob, unwrappedImPermuted, start, end));
                                if ( performTiming_ ) { gt_timer3_.stop(); }
                            }
                        }

                        ro += jobN;
                    }

                    if ( performTiming_ ) { gt_timer3_.start("permuteROTo3rdDimensionFor3DRecon for unwrapped images ... "); }
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::permute3rdDimensionTo1stDimension(unwrappedImPermuted, workOrder3DT->fullkspace_));
                    if ( performTiming_ ) { gt_timer3_.stop(); }
                }
                else
                {
                    for ( n=0; n<(int)N; n++ )
                    {
                        
                    }
                }
            }
            else
            {
                if ( (refN<N) || (refN==1) )
                {
                    hoNDArray<T> kIm(RO, E1, E2, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin());

                    if ( performTiming_ ) { gt_timer3_.start("grappa 3D allocate buffer for unwarpping ... "); }
                    hoNDArrayMemoryManaged<T> buffer3DT_unwrapping(RO, E1, E2, srcCHA, gtPlus_mem_manager_);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( performTiming_ ) { gt_timer3_.start("grappa 3D apply image domain kernel for every channel ... "); }
                    this->applyImageDomainKernelImage(aliasedIm, kIm, buffer3DT_unwrapping, workOrder3DT->fullkspace_);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder3DT->fullkspace_, "unwarppedIm");
                }
                else
                {
                    hoNDArrayMemoryManaged<T> buffer3DT_unwrapping(RO, E1, E2, srcCHA, dstCHA, gtPlus_mem_manager_);

                    hoNDArray<T> complexIm(RO, E1, E2, dstCHA);
                    for ( n=0; n<(int)N; n++ )
                    {
                        hoNDArray<T> kIm(RO, E1, E2, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin()+n*RO*E1*E2*srcCHA*dstCHA);

                        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kIm, "kIm_n");

                        hoNDArray<T> aliasedImN(RO, E1, E2, srcCHA, aliasedIm.begin()+n*RO*E1*E2*srcCHA);

                        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, aliasedImN, "aliasedIm_n");

                        this->applyImageDomainKernelImage(aliasedImN, kIm, buffer3DT_unwrapping, complexIm);
                        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, complexIm, "complexIm_n");

                        memcpy(workOrder3DT->fullkspace_.begin()+n*RO*E1*E2*dstCHA, complexIm.begin(), sizeof(T)*RO*E1*E2*dstCHA);
                    }
                }
            }

            if ( (workOrder3DT->coilMap_->get_size(0)==RO) 
                && (workOrder3DT->coilMap_->get_size(1)==E1) 
                && (workOrder3DT->coilMap_->get_size(2)==E2) 
                && (workOrder3DT->coilMap_->get_size(3)==dstCHA) )
            {
                if ( performTiming_ ) { gt_timer3_.start("grappa 3D coil combination ... "); }
                gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(workOrder3DT->fullkspace_, *workOrder3DT->coilMap_, workOrder3DT->complexIm_);
                if ( performTiming_ ) { gt_timer3_.stop(); }

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder3DT->complexIm_, "combined");
            }

            if ( performTiming_ ) { gt_timer3_.start("grappa 3D go back to kspace ... "); }
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(workOrder3DT->fullkspace_);
            if ( performTiming_ ) { gt_timer3_.stop(); }
        }
        else
        {
            if ( (refN<N) || (refN==1) )
            {
                if ( performTiming_ ) { gt_timer3_.start("grappa 3D test ... "); }
                hoNDArray<T> unmixCoeff(RO, E1, E2, srcCHA, workOrder3DT->unmixingCoeffIm_->begin());
                if ( performTiming_ ) { gt_timer3_.stop(); }

                if ( performTiming_ ) { gt_timer3_.start("grappa 3D apply unmixing coeff ... "); }
                this->applyUnmixCoeffImage(aliasedIm, unmixCoeff, workOrder3DT->complexIm_);
                if ( performTiming_ ) { gt_timer3_.stop(); }

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder3DT->complexIm_, "unwarppedIm");
            }
            else
            {
                for ( n=0; n<(int)N; n++ )
                {
                    hoNDArray<T> unmixCoeff(RO, E1, E2, srcCHA, workOrder3DT->unmixingCoeffIm_->begin()+n*RO*E1*E2*srcCHA);
                    hoNDArray<T> aliasedImN(RO, E1, E2, srcCHA, aliasedIm.begin()+n*RO*E1*E2*srcCHA);
                    hoNDArray<T> unwarppedIm(RO, E1, E2, 1, workOrder3DT->complexIm_.begin()+n*RO*E1*E2);

                    this->applyUnmixCoeffImage(aliasedImN, unmixCoeff, unwarppedIm);

                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, unwarppedIm, "unwarppedIm");
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

        grappa_.gtPlus_mem_manager_ = this->gtPlus_mem_manager_;

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
