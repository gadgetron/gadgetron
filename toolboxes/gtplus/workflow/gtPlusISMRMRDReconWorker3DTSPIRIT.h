/** \file   gtPlusISMRMRDReconWorker3DTSPIRIT.h
    \brief  Implement the 3DT linear SPIRIT reconstruction
    \author Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"
#include "GadgetronTimer.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorker3DT.h"
#include "gtPlusSPIRIT.h"
#include "gtPlusSPIRIT2DTOperator.h"
#include "gtPlusLSQRSolver.h"

#include "GadgetCloudController.h"
#include "GadgetCloudJobMessageReadWrite.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorker3DTSPIRIT : public gtPlusReconWorker3DT<T>
{
public:

    typedef gtPlusReconWorker3DT<T> BaseClass;
    typedef gtPlusReconWorkOrder3DT<T> WorkOrderType;
    typedef typename BaseClass::value_type value_type;

    gtPlusReconWorker3DTSPIRIT() : spirit_kernelIm_permuted_(false), BaseClass() {}
    virtual ~gtPlusReconWorker3DTSPIRIT() {}

    virtual bool performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT);

    virtual bool performCalibPrep(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, WorkOrderType* workOrder3DT);
    virtual bool performCalibImpl(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, WorkOrderType* workOrder3DT, size_t usedN);

    virtual bool performUnwrapping(gtPlusReconWorkOrder3DT<T>* workOrder3DT, const hoNDArray<T>& data);
    virtual bool performUnwarppingImplROPermuted(gtPlusReconWorkOrder<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& ker, hoNDArray<T>& coilMap, hoNDArray<T>& res);
    virtual bool performUnwarppingImpl(gtPlusReconWorkOrder<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res, size_t n);
    virtual bool performUnwarppingImpl(gtPlusReconJob2DT<T>& job);

    virtual bool computeKSpace(gtPlusReconWorkOrder3DT<T>* workOrder3DT);

    virtual bool autoReconParameter(gtPlusReconWorkOrder<T>* workOrder);

    virtual bool splitJob(gtPlusReconWorkOrder3DT<T>* workOrder3DT, size_t& jobN);

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::verbose_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_mem_manager_;

//protected::

    using BaseClass::ref_src_;
    using BaseClass::ref_dst_;
    using BaseClass::data_dst_;
    using BaseClass::ref_coil_map_dst_;
    using BaseClass::startE1_;
    using BaseClass::endE1_;

    gtPlusSPIRIT<T> spirit_;

    bool spirit_kernelIm_permuted_;
};

template <typename T> 
bool gtPlusReconWorker3DTSPIRIT<T>::computeKSpace(gtPlusReconWorkOrder3DT<T>* workOrder3DT)
{
    bool recon_kspace = true;
    return recon_kspace;
}

template <typename T> 
bool gtPlusReconWorker3DTSPIRIT<T>::autoReconParameter(gtPlusReconWorkOrder<T>* workOrder)
{
    gtPlusReconWorkOrder3DT<T>* workOrder3DT = dynamic_cast<gtPlusReconWorkOrder3DT<T>*>(workOrder);
    if ( workOrder3DT == NULL ) return false;

    double acceFactor = workOrder3DT->acceFactorE1_ * workOrder3DT->acceFactorE2_;

    if ( acceFactor>=16 )
    {
        workOrder3DT->spirit_iter_max_ = 150;
        workOrder3DT->spirit_iter_thres_ = 0.0025;
        workOrder3DT->spirit_reg_lamda_ = 0.01;
    }
    else if ( acceFactor>=12 )
    {
        workOrder3DT->spirit_iter_max_ = 100;
        workOrder3DT->spirit_iter_thres_ = 0.0025;
        workOrder3DT->spirit_reg_lamda_ = 0.01;
    }
    else if ( acceFactor>=9 )
    {
        workOrder3DT->spirit_iter_max_ = 100;
        workOrder3DT->spirit_iter_thres_ = 0.0025;
        workOrder3DT->spirit_reg_lamda_ = 0.01;
    }
    else if ( acceFactor>=6 )
    {
        workOrder3DT->spirit_iter_max_ = 100;
        workOrder3DT->spirit_iter_thres_ = 0.0025;
        workOrder3DT->spirit_reg_lamda_ = 0.01;
    }
    else if ( acceFactor>=4 )
    {
        workOrder3DT->spirit_iter_max_ = 70;
        workOrder3DT->spirit_iter_thres_ = 0.005;
        workOrder3DT->spirit_reg_lamda_ = 0.01;
    }
    else
    {
        workOrder3DT->spirit_iter_max_ = 50;
        workOrder3DT->spirit_iter_thres_ = 0.005;
        workOrder3DT->spirit_reg_lamda_ = 0.01;

        if ( workOrder3DT->recon_algorithm_ == ISMRMRD_SPIRIT )
        {
            workOrder3DT->spirit_iter_thres_ = 0.005;
        }
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTSPIRIT<T>::
performCalibPrep(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, WorkOrderType* workOrder3DT)
{
    spirit_.performTiming_ = performTiming_;
    spirit_.debugFolder_ = debugFolder_;

    size_t RO = workOrder3DT->data_.get_size(0);
    size_t E1 = workOrder3DT->data_.get_size(1);
    size_t E2 = workOrder3DT->data_.get_size(2);
    size_t N = workOrder3DT->data_.get_size(4);
    size_t srcCHA = ref_src.get_size(3);

    size_t refRO = ref_dst.get_size(0);
    size_t refE1 = ref_dst.get_size(1);
    size_t refE2 = ref_dst.get_size(2);
    size_t refN = ref_dst.get_size(4);
    size_t dstCHA = ref_dst.get_size(3);

    bool reconKSpace = this->computeKSpace(workOrder3DT);

    size_t kRO = workOrder3DT->spirit_kSize_RO_;
    size_t kE1 = workOrder3DT->spirit_kSize_E1_;
    size_t kE2 = workOrder3DT->spirit_kSize_E2_;

    size_t oRO = workOrder3DT->spirit_oSize_RO_;
    size_t oE1 = workOrder3DT->spirit_oSize_E1_;
    size_t oE2 = workOrder3DT->spirit_oSize_E2_;

    workOrder3DT->kernel_->create(kRO, kE1, kE2, srcCHA, dstCHA, oRO, oE1, oE2, refN);

    size_t jobN;
    bool splitJobs = this->splitJob(workOrder3DT, jobN);

    if ( !splitJobs )
    {
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("allocate image domain kernel ... "));
        if ( gtPlus_mem_manager_ )
        {
            if ( workOrder3DT->kernelIm_->get_number_of_elements() != (size_t)RO*E1*E2*srcCHA*dstCHA*refN )
            {
                workOrder3DT->kernelIm_->create(E1, E2, RO, srcCHA, dstCHA, refN, (T*)(gtPlus_mem_manager_->allocate(sizeof(T)*(size_t)RO*E1*E2*srcCHA*dstCHA*refN)));
            }
        }
        else
        {
            workOrder3DT->kernelIm_->create(E1, E2, RO, srcCHA, dstCHA, refN);
            // pre-set to zero is needed here
            memset(workOrder3DT->kernelIm_->begin(), 0, workOrder3DT->kernelIm_->get_number_of_bytes());
        }
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());
    }
    else
    {
        size_t convKE1 = 2*kE1-1;
        size_t convKE2 = 2*kE2-1;

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("allocate image domain kernel only along RO ... "));
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
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTSPIRIT<T>::
performCalibImpl(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, WorkOrderType* workOrder3DT, size_t usedN)
{
    size_t RO = workOrder3DT->data_.get_size(0);
    size_t E1 = workOrder3DT->data_.get_size(1);
    size_t E2 = workOrder3DT->data_.get_size(2);
    size_t N = workOrder3DT->data_.get_size(4);
    size_t srcCHA = ref_src.get_size(3);

    size_t refRO = ref_dst.get_size(0);
    size_t refE1 = ref_dst.get_size(1);
    size_t refE2 = ref_dst.get_size(2);
    size_t refN = ref_dst.get_size(4);
    size_t dstCHA = ref_dst.get_size(3);

    bool reconKSpace = this->computeKSpace(workOrder3DT);

    size_t kRO = workOrder3DT->spirit_kSize_RO_;
    size_t kE1 = workOrder3DT->spirit_kSize_E1_;
    size_t kE2 = workOrder3DT->spirit_kSize_E2_;

    size_t oRO = workOrder3DT->spirit_oSize_RO_;
    size_t oE1 = workOrder3DT->spirit_oSize_E1_;
    size_t oE2 = workOrder3DT->spirit_oSize_E2_;

    ho4DArray<T> acsSrc(refRO, refE1, refE2, srcCHA, const_cast<T*>(ref_src.begin()+usedN*refRO*refE1*refE2*srcCHA));
    ho4DArray<T> acsDst(refRO, refE1, refE2, dstCHA, const_cast<T*>(ref_dst.begin()+usedN*refRO*refE1*refE2*dstCHA));

    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, acsSrc, "acsSrc");
    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, acsDst, "acsDst");

    hoNDArray<T> ker(kRO, kE1, kE2, srcCHA, dstCHA, oRO, oE1, oE2, workOrder3DT->kernel_->begin()+usedN*kRO*kE1*kE2*srcCHA*dstCHA*oRO*oE1*oE2);

    spirit_.calib_use_gpu_ = workOrder3DT->spirit_use_gpu_;

    GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("SPIRIT 3D calibration ... "));
    GADGET_CHECK_RETURN_FALSE(spirit_.calib3D(acsSrc, acsDst, workOrder3DT->spirit_reg_lamda_, workOrder3DT->spirit_calib_over_determine_ratio_, kRO, kE1, kE2, oRO, oE1, oE2, ker));
    GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, ker, "ker");

    bool minusI = true;

    size_t jobN;
    bool splitJobs = this->splitJob(workOrder3DT, jobN);

    if ( !splitJobs )
    {
        hoNDArray<T> kIm(E1, E2, RO, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin()+usedN*E1*E2*RO*srcCHA*dstCHA);

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("SPIRIT 3D image domain kernel ... "));
        GADGET_CHECK_RETURN_FALSE(spirit_.imageDomainKernel3D(ker, kRO, kE1, kE2, oRO, oE1, oE2, RO, E1, E2, kIm, minusI));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        if ( !debugFolder_.empty() )
        {
            hoNDArray<T> kImACha(E1, E2, RO, srcCHA, kIm.begin());
            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kImACha, "kImACha");
        }
    }
    else
    {
        size_t convKE1 = 2*kE1-1;
        size_t convKE2 = 2*kE2-1;

        hoNDArray<T> kIm(convKE1, convKE2, RO, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin()+usedN*convKE1*convKE2*RO*srcCHA*dstCHA);

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("SPIRIT 3D image domain kernel only along RO ... "));
        GADGET_CHECK_RETURN_FALSE(spirit_.imageDomainKernelRO3D(ker, kRO, kE1, kE2, oRO, oE1, oE2, RO, kIm, minusI));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        if ( !debugFolder_.empty() )
        {
            hoNDArray<T> kImROACha(convKE1, convKE2, RO, srcCHA, kIm.begin());
            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kImROACha, "kImROACha");
        }
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTSPIRIT<T>::
splitJob(gtPlusReconWorkOrder3DT<T>* workOrder3DT, size_t& jobN)
{
    size_t RO = workOrder3DT->data_.get_size(0);
    size_t E1 = workOrder3DT->data_.get_size(1);
    size_t E2 = workOrder3DT->data_.get_size(2);

    size_t srcCHA = workOrder3DT->kernelIm_->get_size(3);
    size_t dstCHA = workOrder3DT->kernelIm_->get_size(4);

    bool splitByS = workOrder3DT->job_split_by_S_;
    jobN = workOrder3DT->job_num_of_N_;
    size_t jobMegaBytes = workOrder3DT->job_max_Megabytes_;

    bool splitJobs = (splitByS==true || jobN>0);
    if ( !splitJobs )
    {
        if ( jobMegaBytes>0 )
        {
            size_t jobN = jobMegaBytes/(E1*E2*srcCHA*dstCHA*sizeof(T)/1024/1024);
            if ( jobN < RO ) splitJobs = true;
            GADGET_MSG("SPIRIT - 3DT - size of largest job : " << jobN);
        }
    }
    if ( jobN >= RO ) splitJobs = false;

    return splitJobs;
}

template <typename T> 
bool gtPlusReconWorker3DTSPIRIT<T>::
performUnwrapping(gtPlusReconWorkOrder3DT<T>* workOrder3DT, const hoNDArray<T>& data_dst)
{
    try
    {
        int n;

        size_t RO = workOrder3DT->data_.get_size(0);
        size_t E1 = workOrder3DT->data_.get_size(1);
        size_t E2 = workOrder3DT->data_.get_size(2);
        size_t N = workOrder3DT->data_.get_size(4);

        size_t kImE1 = workOrder3DT->kernelIm_->get_size(0);
        size_t kImE2 = workOrder3DT->kernelIm_->get_size(1);
        size_t kImRO = workOrder3DT->kernelIm_->get_size(2);
        size_t srcCHA = workOrder3DT->kernelIm_->get_size(3);
        size_t dstCHA = workOrder3DT->kernelIm_->get_size(4);

        size_t refN = workOrder3DT->kernelIm_->get_size(5);

        workOrder3DT->complexIm_.create(RO, E1, E2, N);

        // downstream coil compression is not supported here
        // kspace is always reconed
        workOrder3DT->fullkspace_ = data_dst;

        // compute the scaling factor
        typename realType<T>::Type scaleFactor = 1.0;
        hoNDArray<T> kspaceForScaleFactor(RO, E1, E2, srcCHA, const_cast<T*>(data_dst.begin()));
        Gadgetron::norm2(kspaceForScaleFactor, scaleFactor);
        scaleFactor /= (value_type)(RO*std::sqrt(double(srcCHA)));

        workOrder3DT->spirit_ncg_scale_factor_ = scaleFactor;

        size_t indMax;
        hoNDArray<value_type> mag;
        Gadgetron::abs(kspaceForScaleFactor, mag);
        value_type maxMag;
        Gadgetron::maxAbsolute(mag, maxMag, indMax);
        workOrder3DT->spirit_slep_scale_factor_ = maxMag;

        // split the jobs
        size_t jobMegaBytes = workOrder3DT->job_max_Megabytes_;
        size_t jobN = workOrder3DT->job_num_of_N_;
        bool splitJobs = this->splitJob(workOrder3DT, jobN);
        size_t maxNumOfBytesPerJob = jobMegaBytes*1024*1024;

        size_t overlapN = workOrder3DT->job_overlap_;
        if ( workOrder3DT->recon_algorithm_==ISMRMRD_SPIRIT )
        {
            overlapN = 0;
        }

        if ( splitJobs )
        {
            // hoNDArrayMemoryManaged<T> kspaceIfftRO(RO, E1, E2, srcCHA, N, gtPlus_mem_manager_);
            hoNDArray<T> kspaceIfftRO(RO, E1, E2, srcCHA, N);
            GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(data_dst, kspaceIfftRO));
            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kspaceIfftRO, "kspaceIfftRO");

            // hoNDArrayMemoryManaged<T> kspaceIfftROPermuted(E1, E2, srcCHA, RO, N, gtPlus_mem_manager_);
            hoNDArray<T> kspaceIfftROPermuted(E1, E2, srcCHA, RO, N);
            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("permute kspace RO to 4th dimension ... "));
            GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo4thDimensionFor3DRecon(kspaceIfftRO, kspaceIfftROPermuted));
            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());
            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kspaceIfftROPermuted, "kspaceIfftROPermuted");

            hoNDArrayMemoryManaged<T> kerPermuted;
            if ( !spirit_kernelIm_permuted_ )
            {
                spirit_kernelIm_permuted_ = true;

                size_t kerN = kImE1*kImE2*srcCHA*dstCHA*kImRO*N;
                size_t kerImSize = sizeof(T)*kerN;
                GADGET_MSG("SPIRIT - 3DT - image domain kernel size : " << kerImSize/1024.0/1024 << " MBytes ... ");
                size_t maxFreeChunk = gtPlus_mem_manager_->maxFreeMemoryChunkSize();
                GADGET_MSG("SPIRIT - 3DT - maximal free chunk of managed memory : " << maxFreeChunk/1024.0/1024 << " MBytes ... ");

                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("allocate permuted kernel ... "));
                if ( maxFreeChunk >= kerImSize )
                {
                    kerPermuted.setMemoryManager(gtPlus_mem_manager_);
                    kerPermuted.create(kImE1, kImE2, srcCHA, dstCHA, kImRO, N);
                }
                else
                {
                    GADGET_MSG("use unmanaged memory ... ");
                    T* pData = new T[kerN];
                    kerPermuted.create(kImE1, kImE2, srcCHA, dstCHA, kImRO, N, pData, true);
                }
                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *workOrder3DT->kernelIm_, "kernelImBeforePermuted");

                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("permute kernel RO to 5th dimension ... "));
                // GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo5thDimensionFor3DRecon( *workOrder3DT->kernelIm_, kerPermuted));
                GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteE2To5thDimension( *workOrder3DT->kernelIm_, kerPermuted));
                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kerPermuted, "kerPermuted");

                workOrder3DT->kernelIm_->reshape(kerPermuted.get_dimensions());
                *workOrder3DT->kernelIm_ = kerPermuted;

                kerPermuted.clear();

                kerPermuted.create(kImE1, kImE2, srcCHA, dstCHA, kImRO, N, workOrder3DT->kernelIm_->begin());
            }
            else
            {
                kerPermuted.create(E1, E2, srcCHA, dstCHA, RO, N, workOrder3DT->kernelIm_->begin());
            }

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kerPermuted, "kerPermuted_Used");

            gtPlusReconWorkOrder3DT<T> workOrder3DTJobSplit;
            workOrder3DT->duplicate(workOrder3DTJobSplit);

            boost::shared_ptr< hoNDArray<T> > coilMapPermuted = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>()) ;
            if ( workOrder3DT->coilMap_->get_number_of_elements() > 0 )
            {
                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("permute coil map RO to 4th dimension ... "));
                GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo4thDimensionFor3DRecon(*workOrder3DT->coilMap_, *coilMapPermuted));
                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());
                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kspaceIfftROPermuted, "coilMapPermuted");

                workOrder3DTJobSplit.coilMap_ = coilMapPermuted;
            }

            bool runJobsOnCloud = workOrder3DT->CloudComputing_;
            unsigned int cloudSize = workOrder3DT->CloudSize_;
            bool runJobsOnLocalNode = workOrder3DT->job_perform_on_control_node_;

            std::vector<gtPlusReconJob2DT<T> > jobList;

            if ( runJobsOnCloud )
            {
                unsigned int j;

                GADGET_CHECK_RETURN_FALSE(this->estimateJobSize(workOrder3DT, maxNumOfBytesPerJob, overlapN, cloudSize, jobN));

                //GADGET_MSG("SPIRIT - 3DT - cloudSize is " << cloudSize << " - RO is " << RO << " ... ");
                //unsigned int nodeN = cloudSize;
                //if ( runJobsOnLocalNode ) nodeN++;
                //GADGET_MSG("SPIRIT - 3DT - runJobsOnLocalNode is " << runJobsOnLocalNode << " - nodeN is " << nodeN << " - overlapN is " << overlapN << " ... ");

                //// adjust jobN according to cloud size
                //jobN = std::ceil( (double)(RO+overlapN*(nodeN-1))/(double)nodeN );

                //size_t numOfBytesPerJob = sizeof(T)*( E1*E2*srcCHA*dstCHA*jobN + 2*E1*E2*srcCHA*jobN );

                //while ( numOfBytesPerJob > 2.2*1024*1024*1024-64.0*1024*1024 )
                //{
                //    nodeN *= 2;
                //    jobN = std::ceil( (double)(RO+overlapN*(nodeN-1))/(double)nodeN );
                //    numOfBytesPerJob = sizeof(T)*( E1*E2*srcCHA*dstCHA*jobN + 2*E1*E2*srcCHA*jobN );
                //}

                //GADGET_MSG("SPIRIT - 3DT - every job will have " << numOfBytesPerJob/1024.0/1024 << " MBytes ... ");

                // split the job
                GADGET_CHECK_RETURN_FALSE(this->splitReconJob(&workOrder3DTJobSplit, kspaceIfftROPermuted, kerPermuted, workOrder3DT->job_split_by_S_, jobN, jobMegaBytes, overlapN, jobList));

                std::vector<gtPlusReconJob2DT<T> > completedJobList(jobList.size());

                for ( j=0; j<jobList.size(); j++ )
                {
                    jobList[j].workOrder2DT.duplicate(completedJobList[j].workOrder2DT);
                    completedJobList[j].job_index_startN_ = jobList[j].job_index_startN_;
                    completedJobList[j].job_index_endN_ = jobList[j].job_index_endN_;
                    completedJobList[j].job_index_S_ = jobList[j].job_index_S_;
                }

                GADGET_MSG("SPIRIT - 3DT - total job : " << jobList.size() << " - job N : " << jobN << " - cloud size : " << cloudSize);

                unsigned int numOfJobRunOnCloud = (unsigned int)(jobList.size() - jobList.size()/(cloudSize+1));
                if ( !runJobsOnLocalNode ) numOfJobRunOnCloud = (unsigned int)jobList.size();

                typedef Gadgetron::GadgetCloudController< gtPlusReconJob2DT<T> > GTCloudControllerType;
                GTCloudControllerType controller;

                if (controller.open () == -1)
                {
                    GADGET_ERROR_MSG("Cloud controller cannot open the cloud ...");
                    controller.handle_close (ACE_INVALID_HANDLE, 0);
                    runJobsOnCloud = false;
                }
                else
                {
                    std::vector<gtPlusReconJob2DT<T>* > jobListCloud(numOfJobRunOnCloud);
                    std::vector<gtPlusReconJob2DT<T>* > completedJobListCloud(numOfJobRunOnCloud);
                    std::vector<int> node_ids(numOfJobRunOnCloud);

                    GADGET_CHECK_RETURN_FALSE(this->scheduleJobForNodes(workOrder3DT, numOfJobRunOnCloud, node_ids));

                    for ( j=0; j<numOfJobRunOnCloud; j++ )
                    {
                        // node_ids[j] = j%cloudSize;
                        jobListCloud[j] = &jobList[j];
                        completedJobListCloud[j] = &completedJobList[j];
                        GADGET_MSG("--> job " << j << " runs on node " << node_ids[j] << " ... ");
                    }

                    std::vector<GadgetMessageReader*> readers(cloudSize, NULL);
                    std::vector<GadgetMessageWriter*> writers(cloudSize, NULL);

                    for ( j=0; j<cloudSize; j++ )
                    {
                        readers[j] = new GtPlusCloudJobMessageReaderCPFL();
                        writers[j] = new GtPlusCloudJobMessageWriterCPFL();
                    }

                    if ( controller.createConnector(workOrder3DT->gt_cloud_, GADGET_MESSAGE_CLOUD_JOB, readers, GADGET_MESSAGE_CLOUD_JOB, writers) != 0 )
                    {
                        GADGET_ERROR_MSG("Cloud controller creates connectors failed ...");
                        controller.handle_close (ACE_INVALID_HANDLE, 0);
                        runJobsOnCloud = false;
                    }
                    else if ( controller.connectToCloud(workOrder3DT->gt_cloud_) != 0 )
                    {
                        GADGET_ERROR_MSG("Cloud controller cannot connect to the cloud ...");
                        controller.handle_close (ACE_INVALID_HANDLE, 0);
                        runJobsOnCloud = false;
                    }
                    else
                    {
                        if ( controller.runJobsOnCloud(jobListCloud, completedJobListCloud, node_ids) != 0 )
                        {
                            GADGET_ERROR_MSG("Cloud controller runs jobs on the cloud failed ...");
                            controller.closeCloudNode();
                            controller.handle_close (ACE_INVALID_HANDLE, 0);
                            runJobsOnCloud = false;
                        }
                        else
                        {
                            controller.closeCloudNode();

                            // run the left over jobs on the local computer
                            for ( j=numOfJobRunOnCloud; j<jobList.size(); j++ )
                            {
                                GADGET_MSG("SPIRIT - 3DT - job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("SPIRIT 3DT ... "));
                                GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(jobList[j]));
                                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

                                std::ostringstream ostr;
                                ostr << "job_fullkspace" << "_" << j;
                                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, jobList[j].res, ostr.str());
                            }

                            // wait the cloud job to complete
                            controller.waitForJobToComplete();

                            // combine results from cloud and local run
                            for ( j=0; j<numOfJobRunOnCloud; j++ )
                            {
                                jobList[j].res = controller.completed_job_list_[j]->res;
                                jobList[j].complexIm = controller.completed_job_list_[j]->complexIm;
                            }

                            // if some jobs are not actually completed, process them
                            for ( j=0; j<numOfJobRunOnCloud; j++ )
                            {
                                if ( 
                                    !jobList[j].res.dimensions_equal(&jobList[j].kspace) 
                                        && 
                                    ( jobList[j].complexIm.get_size(0)!= jobList[j].kspace.get_size(0) 
                                    || jobList[j].complexIm.get_size(1)!= jobList[j].kspace.get_size(1) 
                                    || jobList[j].complexIm.get_size(2)!= jobList[j].kspace.get_size(2) ) 
                                   )
                                {
                                    GADGET_MSG("SPIRIT - 3DT - uncompleted cloud job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                                    GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("SPIRIT 3DT ... "));
                                    GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(jobList[j]));
                                    GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

                                    std::ostringstream ostr;
                                    ostr << "job_fullkspace" << "_" << j;
                                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, jobList[j].res, ostr.str());
                                }
                            }
                        }
                    }
                }
            }

            if ( !runJobsOnCloud )
            {
                // split the job
                GADGET_CHECK_RETURN_FALSE(this->splitReconJob(&workOrder3DTJobSplit, kspaceIfftROPermuted, kerPermuted, workOrder3DT->job_split_by_S_, jobN, jobMegaBytes, overlapN, jobList));

                GADGET_MSG("SPIRIT - 3DT - total job : " << jobList.size());

                size_t j;
                for ( j=0; j<jobList.size(); j++ )
                {
                    GADGET_MSG("SPIRIT - 3DT - job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                    GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("SPIRIT 3DT ... "));
                    GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(jobList[j]));
                    GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

                    std::ostringstream ostr;
                    ostr << "job_fullkspace" << "_" << j;
                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, jobList[j].res, ostr.str());
                }
            }

            // combine the job
            workOrder3DTJobSplit.fullkspace_.create(E1, E2, dstCHA, RO, N);
            GADGET_CHECK_RETURN_FALSE(this->combineReconJob(&workOrder3DTJobSplit, jobList, RO, N));

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder3DTJobSplit.fullkspace_, "job_combined_fullkspace");

            // clear the memory
            jobList.clear();

            // permute the unwrapped kspace
            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("permtue RO to 1st dimension ... "));
            GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo1stDimensionFor3DRecon(workOrder3DTJobSplit.fullkspace_, kspaceIfftRO));
            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kspaceIfftRO, "res_fullkspace_ROinIm");

            // perform fft along the first dimension
            GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(kspaceIfftRO, workOrder3DT->fullkspace_));
            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder3DT->fullkspace_, "res_3DSpirit");
        }
        else
        {
            for ( n=0; n<(int)N; n++ )
            {
                size_t kernelN = n;
                if ( kernelN >= refN ) kernelN = refN-1;

                hoNDArray<T> kIm(E1, E2, RO, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin()+kernelN*RO*E1*E2*srcCHA*dstCHA);

                hoNDArray<T> aliasedKSpace(RO, E1, E2, srcCHA, const_cast<T*>(data_dst.begin())+n*RO*E1*E2*srcCHA);

                hoNDArray<T> unwarppedKSpace(RO, E1, E2, dstCHA, workOrder3DT->fullkspace_.begin()+n*RO*E1*E2*dstCHA);

                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D unwrapping ... "));
                GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(workOrder3DT, aliasedKSpace, kIm, unwarppedKSpace, n));
                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, unwarppedKSpace, "unwarppedKSpace");
            }
        }

        hoNDArrayMemoryManaged<T> complexImMultiChannel(RO, E1, E2, dstCHA, N, gtPlus_mem_manager_);

        if ( (workOrder3DT->coilMap_->get_size(0)==RO) 
            && (workOrder3DT->coilMap_->get_size(1)==E1) 
            && (workOrder3DT->coilMap_->get_size(2)==E2) 
            && (workOrder3DT->coilMap_->get_size(3)==dstCHA) )
        {
            hoNDArrayMemoryManaged<T> complexImMultiChannel(RO, E1, E2, dstCHA, N, gtPlus_mem_manager_);
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->fullkspace_, complexImMultiChannel);
            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, complexImMultiChannel, "unwarppedComplexIm");

            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D coil combination ... "));
            gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(complexImMultiChannel, *workOrder3DT->coilMap_, workOrder3DT->complexIm_);
            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder3DT->complexIm_, "combined");
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker3DTSPIRIT<T>::performUnwrapping(gtPlusReconWorkOrder3DT<T>* workOrder3DT, const hoNDArray<T>& data) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTSPIRIT<T>::
performUnwarppingImplROPermuted(gtPlusReconWorkOrder<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& ker, hoNDArray<T>& /*coilMap*/, hoNDArray<T>& res)
{
    try
    {
        size_t E1 = kspace.get_size(0);
        size_t E2 = kspace.get_size(1);
        size_t RO = kspace.get_size(3);

        size_t kerE1 = ker.get_size(0);
        size_t kerE2 = ker.get_size(1);
        size_t srcCHA = ker.get_size(2);
        size_t dstCHA = ker.get_size(3);
        size_t kerN = ker.get_size(5);

        hoNDArray<T>* kerIm = &ker;
        hoNDArray<T> kerImE1E2RO;
        if ( kerE1!=E1 || kerE2!=E2 )
        {
            GADGET_MSG("gtPlusReconWorker3DTSPIRIT, kerE1!=E1 || kerE2!=E2, kernel needs to be converted along E1 and E2 ... ");

            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("kernel conversion along E1 and E2 ... "));

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

            GADGET_CHECK_RETURN_FALSE(spirit_.imageDomainKernelE1E2RO(ker, (int)E1, (int)E2, kerImE1E2RO));
            kerIm = &kerImE1E2RO;

            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());
        }

        res.create(kspace.get_dimensions());

        long long NUM = (long long)RO;

        #ifdef USE_OMP
            int numThreads = NUM;

            int numOpenMPProcs = omp_get_num_procs();
            GADGET_MSG("gtPlusReconWorker3DTSPIRIT, numOpenMPProcs : " << numOpenMPProcs);

            if ( numThreads > numOpenMPProcs ) numThreads = numOpenMPProcs;
            GADGET_MSG("gtPlusReconWorker3DTSPIRIT, numThreads : " << numThreads);
        #endif

        long long t;

        hoNDArray<T> ker_Shifted(kerIm);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(*kerIm, ker_Shifted);
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, ker_Shifted, "ker_Shifted");

        hoNDArray<T> kspace_Shifted(kspace);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(kspace, kspace_Shifted);
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kspace_Shifted, "kspace_Shifted");

        #ifdef GCC_OLD_FLAG
            #pragma omp parallel default(none) private(t) shared(RO, E1, E2, srcCHA, dstCHA, workOrder3DT, kspace_Shifted, ker_Shifted, NUM) if ( NUM > 1 ) num_threads( numThreads )
        #else
            #pragma omp parallel default(none) private(t) shared(RO, E1, E2, srcCHA, dstCHA, workOrder3DT, NUM, kspace_Shifted, ker_Shifted, res) if ( NUM > 1 ) num_threads( numThreads )
        #endif
        {
            gtPlusSPIRIT2DOperator<T> spirit;
            spirit.setMemoryManager(gtPlus_mem_manager_);
            spirit.use_symmetric_spirit_ = false;
            spirit.use_non_centered_fft_ = true;

            hoNDArray<T> x0(E1, E2, srcCHA);
            Gadgetron::clear(x0);

            gtPlusLinearSolver<hoNDArray<T>, hoNDArray<T>, gtPlusSPIRIT2DOperator<T> >* pCGSolver;

            pCGSolver = new gtPlusLSQRSolver<hoNDArray<T>, hoNDArray<T>, gtPlusSPIRIT2DOperator<T> >();

            gtPlusLinearSolver<hoNDArray<T>, hoNDArray<T>, gtPlusSPIRIT2DOperator<T> >& cgSolver = *pCGSolver;

            cgSolver.iterMax_ = workOrder3DT->spirit_iter_max_;
            cgSolver.thres_ = (value_type)workOrder3DT->spirit_iter_thres_;
            cgSolver.printIter_ = workOrder3DT->spirit_print_iter_;

            cgSolver.set(spirit);

            hoNDArray<T> b(E1, E2, srcCHA);

            #pragma omp for
            for ( t=0; t<NUM; t++ )
            {
                size_t ro = t;

                hoNDArray<T> kspaceCurr(E1, E2, srcCHA, kspace_Shifted.begin()+ro*E1*E2*srcCHA);
                hoNDArray<T> resCurr(E1, E2, dstCHA, res.begin()+ro*E1*E2*dstCHA);

                // solve the 2D spirit problem
                Gadgetron::clear(x0);

                boost::shared_ptr<hoNDArray<T> > kerCurr(new hoNDArray<T>(E1, E2, srcCHA, dstCHA, ker_Shifted.begin()+ro*E1*E2*srcCHA*dstCHA));

                spirit.setForwardKernel(kerCurr, false);

                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(E1, E2, srcCHA, kspaceCurr.begin()));
                spirit.setAcquiredPoints(acq);

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *kerCurr, "spirit3D_ker");
                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *acq, "spirit3D_kspace");

                cgSolver.x0_ = acq.get();

                // compute rhs
                spirit.computeRighHandSide(*acq, b);

                // solve
                cgSolver.solve(b, resCurr);

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, resCurr, "unwarppedKSpace_t");

                // restore the acquired points
                spirit.restoreAcquiredKSpace(*acq, resCurr);

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, resCurr, "unwarppedKSpace_t_setAcq");
            }

            delete pCGSolver;
        }

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res, "res_Shifted");

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fftshift2D(res, kspace_Shifted);
        res = kspace_Shifted;

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res, "resPermuted");
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker3DTSPIRIT<T>::performUnwarppingImplROPermuted(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTSPIRIT<T>::
performUnwarppingImpl(gtPlusReconWorkOrder<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res, size_t n)
{
    try
    {
        // RO, E1, E2, srcCHA, dstCHA, N
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t E2 = kspace.get_size(2);

        size_t srcCHA = adj_forward_G_I.get_size(3);
        size_t dstCHA = adj_forward_G_I.get_size(4);

        // perform the 3D recon by read-out decoupling
        hoNDArrayMemoryManaged<T> resDecoupled(E1, E2, dstCHA, RO, gtPlus_mem_manager_);

        hoNDArrayMemoryManaged<T> kspaceIfftRO(RO, E1, E2, srcCHA, gtPlus_mem_manager_);
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(kspace, kspaceIfftRO));
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kspaceIfftRO, "kspaceIfftRO");

        hoNDArrayMemoryManaged<T> kspaceIfftROPermuted(E1, E2, srcCHA, RO, gtPlus_mem_manager_);

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("permtue RO to 4th dimension ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo4thDimensionFor3DRecon(kspaceIfftRO, kspaceIfftROPermuted));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kspaceIfftROPermuted, "kspaceIfftROPermuted");

        T* pKspaceIfftROPermuted = kspaceIfftROPermuted.begin();

        T* pG_I = adj_forward_G_I.begin();

        long long NUM = (long long)RO;

        long long t;

        #pragma omp parallel default(none) private(t) shared(RO, E1, E2, srcCHA, dstCHA, workOrder3DT, NUM, resDecoupled, pKspaceIfftROPermuted, pG_I) if ( NUM > 6 ) num_threads( (int)((NUM<16) ? NUM : 16) )
        {
            hoNDArrayMemoryManaged<T> adjForG_I_Decoupled(E1, E2, srcCHA, dstCHA, gtPlus_mem_manager_);
            T* pDecoupledG_I = adjForG_I_Decoupled.begin();

            gtPlusSPIRIT2DOperator<T> spirit;
            spirit.setMemoryManager(gtPlus_mem_manager_);
            spirit.use_symmetric_spirit_ = false;

            hoNDArray<T> x0(E1, E2, srcCHA);
            Gadgetron::clear(x0);

            gtPlusLinearSolver<hoNDArray<T>, hoNDArray<T>, gtPlusSPIRIT2DOperator<T> >* pCGSolver;
            pCGSolver = new gtPlusLSQRSolver<hoNDArray<T>, hoNDArray<T>, gtPlusSPIRIT2DOperator<T> >();
            gtPlusLinearSolver<hoNDArray<T>, hoNDArray<T>, gtPlusSPIRIT2DOperator<T> >& cgSolver = *pCGSolver;

            cgSolver.iterMax_ = workOrder3DT->spirit_iter_max_;
            cgSolver.thres_ = (value_type)workOrder3DT->spirit_iter_thres_;
            cgSolver.printIter_ = workOrder3DT->spirit_print_iter_;

            cgSolver.set(spirit);

            hoNDArray<T> b(E1, E2, srcCHA);

            #pragma omp for
            for ( t=0; t<NUM; t++ )
            {
                size_t ro = t;

                hoNDArray<T> kspace_DeDecoupled(E1, E2, srcCHA, pKspaceIfftROPermuted+ro*E1*E2*srcCHA);
                hoNDArray<T> resCurr(E1, E2, dstCHA, resDecoupled.begin()+ro*E1*E2*dstCHA);

                // fill in kernel and kspace
                size_t scha, dcha;

                for ( dcha=0; dcha<dstCHA; dcha++)
                {
                    for ( scha=0; scha<srcCHA; scha++)
                    {

                        T* pDst = pDecoupledG_I + scha*E1*E2+dcha*E1*E2*srcCHA;
                        T* pSrc = pG_I + ro*E1*E2+scha*RO*E1*E2+dcha*RO*E1*E2*srcCHA;
                        memcpy(pDst, pSrc, sizeof(T)*E1*E2);
                    }
                }

                // solve the 2D spirit problem
                Gadgetron::clear(x0);

                boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(E1, E2, srcCHA, dstCHA, pDecoupledG_I));

                spirit.setForwardKernel(ker, false);

                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(E1, E2, srcCHA, kspace_DeDecoupled.begin()));
                spirit.setAcquiredPoints(acq);

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *ker, "spirit3D_ker");
                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *acq, "spirit3D_kspace");

                cgSolver.x0_ = acq.get();

                // compute rhs
                spirit.computeRighHandSide(*acq, b);

                // solve
                cgSolver.solve(b, resCurr);

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, resCurr, "unwarppedKSpace_t");

                // restore the acquired points
                spirit.restoreAcquiredKSpace(*acq, resCurr);

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, resCurr, "unwarppedKSpace_t_setAcq");
            }

            delete pCGSolver;
        }

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, resDecoupled, "resDecoupled");

        // permute the unwrapped kspace
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("permtue RO to 1st dimension ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo1stDimensionFor3DRecon(resDecoupled, kspaceIfftRO));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        // perform fft along the first dimension
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(kspaceIfftRO, res));
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res, "res_3DSpirit");
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker3DTSPIRIT<T>::performUnwarppingImpl(gtPlusReconWorkOrder3DT<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTSPIRIT<T>::
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
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker3DTSPIRIT<T>::performUnwarppingImpl(job) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DTSPIRIT<T>::performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(workOrder3DT!=NULL);

        spirit_.gtPlus_mem_manager_ = this->gtPlus_mem_manager_;

        // call the BaseClass
        GADGET_CHECK_RETURN_FALSE(BaseClass::performRecon(workOrder3DT));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker3DTSPIRIT<T>::performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT) ... ");
        return false;
    }

    return true;
}

}}
