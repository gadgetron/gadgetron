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
#include "hoLSQRSolver.h"
#include "hoSPIRIT2DOperator.h"

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
        if ( performTiming_ ) { gt_timer3_.start("allocate image domain kernel ... "); }
        workOrder3DT->kernelIm_->create(E1, E2, RO, srcCHA, dstCHA, refN);
        // pre-set to zero is needed here
        memset(workOrder3DT->kernelIm_->begin(), 0, workOrder3DT->kernelIm_->get_number_of_bytes());
        if ( performTiming_ ) { gt_timer3_.stop(); }
    }
    else
    {
        size_t convKE1 = 2*kE1-1;
        size_t convKE2 = 2*kE2-1;

        if ( performTiming_ ) { gt_timer3_.start("allocate image domain kernel only along RO ... "); }
        workOrder3DT->kernelIm_->create(convKE1, convKE2, RO, srcCHA, dstCHA, refN);
        // pre-set to zero is needed here
        memset(workOrder3DT->kernelIm_->begin(), 0, workOrder3DT->kernelIm_->get_number_of_bytes());
        if ( performTiming_ ) { gt_timer3_.stop(); }
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

    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(acsSrc, debugFolder_+"acsSrc"); }
    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(acsDst, debugFolder_+"acsDst"); }

    hoNDArray<T> ker(kRO, kE1, kE2, srcCHA, dstCHA, oRO, oE1, oE2, workOrder3DT->kernel_->begin()+usedN*kRO*kE1*kE2*srcCHA*dstCHA*oRO*oE1*oE2);

    spirit_.calib_use_gpu_ = workOrder3DT->spirit_use_gpu_;

    if ( performTiming_ ) { gt_timer3_.start("SPIRIT 3D calibration ... "); }
    GADGET_CHECK_RETURN_FALSE(spirit_.calib3D(acsSrc, acsDst, workOrder3DT->spirit_reg_lamda_, workOrder3DT->spirit_calib_over_determine_ratio_, kRO, kE1, kE2, oRO, oE1, oE2, ker));
    if ( performTiming_ ) { gt_timer3_.stop(); }

    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ker, debugFolder_+"ker"); }

    bool minusI = true;

    size_t jobN;
    bool splitJobs = this->splitJob(workOrder3DT, jobN);

    if ( !splitJobs )
    {
        hoNDArray<T> kIm(E1, E2, RO, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin()+usedN*E1*E2*RO*srcCHA*dstCHA);

        if ( performTiming_ ) { gt_timer3_.start("SPIRIT 3D image domain kernel ... "); }
        GADGET_CHECK_RETURN_FALSE(spirit_.imageDomainKernel3D(ker, kRO, kE1, kE2, oRO, oE1, oE2, RO, E1, E2, kIm, minusI));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() )
        {
            hoNDArray<T> kImACha(E1, E2, RO, srcCHA, kIm.begin());
            gt_exporter_.exportArrayComplex(kImACha, debugFolder_+"kImACha");
        }
    }
    else
    {
        size_t convKE1 = 2*kE1-1;
        size_t convKE2 = 2*kE2-1;

        hoNDArray<T> kIm(convKE1, convKE2, RO, srcCHA, dstCHA, workOrder3DT->kernelIm_->begin()+usedN*convKE1*convKE2*RO*srcCHA*dstCHA);

        if ( performTiming_ ) { gt_timer3_.start("SPIRIT 3D image domain kernel only along RO ... "); }
        GADGET_CHECK_RETURN_FALSE(spirit_.imageDomainKernelRO3D(ker, kRO, kE1, kE2, oRO, oE1, oE2, RO, kIm, minusI));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() )
        {
            hoNDArray<T> kImROACha(convKE1, convKE2, RO, srcCHA, kIm.begin());
            gt_exporter_.exportArrayComplex(kImROACha, debugFolder_+"kImROACha");
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
            GDEBUG_STREAM("SPIRIT - 3DT - size of largest job : " << jobN);
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
            hoNDArray<T> kspaceIfftRO(RO, E1, E2, srcCHA, N);
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(data_dst, kspaceIfftRO);
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIfftRO, debugFolder_+"kspaceIfftRO"); }

            hoNDArray<T> kspaceIfftROPermuted(E1, E2, srcCHA, RO, N);
            if ( performTiming_ ) { gt_timer3_.start("permute kspace RO to 4th dimension ... "); }

            std::vector<size_t> dim_order(4);
            dim_order[0] = 1;
            dim_order[1] = 2;
            dim_order[2] = 3;
            dim_order[3] = 0;

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&kspaceIfftRO, &kspaceIfftROPermuted, &dim_order));

            if ( performTiming_ ) { gt_timer3_.stop(); }
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIfftROPermuted, debugFolder_+"kspaceIfftROPermuted"); }

            hoNDArray<T> kerPermuted;
            if ( !spirit_kernelIm_permuted_ )
            {
                spirit_kernelIm_permuted_ = true;

                size_t kerN = kImE1*kImE2*srcCHA*dstCHA*kImRO*N;
                size_t kerImSize = sizeof(T)*kerN;
                GDEBUG_STREAM("SPIRIT - 3DT - image domain kernel size : " << kerImSize/1024.0/1024 << " MBytes ... ");

                if ( performTiming_ ) { gt_timer3_.start("allocate permuted kernel ... "); }
                kerPermuted.create(kImE1, kImE2, srcCHA, dstCHA, kImRO, N);
                if ( performTiming_ ) { gt_timer3_.stop(); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder3DT->kernelIm_, debugFolder_+"kernelImBeforePermuted"); }

                if ( performTiming_ ) { gt_timer3_.start("permute kernel RO to 5th dimension ... "); }

                std::vector<size_t> dim_order(5);
                dim_order[0] = 0;
                dim_order[1] = 1;
                dim_order[2] = 3;
                dim_order[3] = 4;
                dim_order[4] = 2;

                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(workOrder3DT->kernelIm_.get(), &kerPermuted, &dim_order));

                if ( performTiming_ ) { gt_timer3_.stop(); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kerPermuted, debugFolder_+"kerPermuted"); }

                workOrder3DT->kernelIm_->reshape(kerPermuted.get_dimensions());
                *workOrder3DT->kernelIm_ = kerPermuted;

                kerPermuted.clear();

                kerPermuted.create(kImE1, kImE2, srcCHA, dstCHA, kImRO, N, workOrder3DT->kernelIm_->begin());
            }
            else
            {
                kerPermuted.create(E1, E2, srcCHA, dstCHA, RO, N, workOrder3DT->kernelIm_->begin());
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kerPermuted, debugFolder_+"kerPermuted_Used"); }

            gtPlusReconWorkOrder3DT<T> workOrder3DTJobSplit;
            workOrder3DT->duplicate(workOrder3DTJobSplit);

            boost::shared_ptr< hoNDArray<T> > coilMapPermuted = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>()) ;
            if ( workOrder3DT->coilMap_->get_number_of_elements() > 0 )
            {
                if ( performTiming_ ) { gt_timer3_.start("permute coil map RO to 4th dimension ... "); }

                coilMapPermuted->create(workOrder3DT->coilMap_->get_size(1), 
                                        workOrder3DT->coilMap_->get_size(2), 
                                        workOrder3DT->coilMap_->get_size(3), 
                                        workOrder3DT->coilMap_->get_size(0), 
                                        workOrder3DT->coilMap_->get_size(4) );

                std::vector<size_t> dim_order(4);
                dim_order[0] = 1;
                dim_order[1] = 2;
                dim_order[2] = 3;
                dim_order[3] = 0;

                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(workOrder3DT->coilMap_.get(), coilMapPermuted.get(), &dim_order));

                if ( performTiming_ ) { gt_timer3_.stop(); }
                if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(*coilMapPermuted, debugFolder_ + "coilMapPermuted"); }

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

                //GDEBUG_STREAM("SPIRIT - 3DT - cloudSize is " << cloudSize << " - RO is " << RO << " ... ");
                //unsigned int nodeN = cloudSize;
                //if ( runJobsOnLocalNode ) nodeN++;
                //GDEBUG_STREAM("SPIRIT - 3DT - runJobsOnLocalNode is " << runJobsOnLocalNode << " - nodeN is " << nodeN << " - overlapN is " << overlapN << " ... ");

                //// adjust jobN according to cloud size
                //jobN = std::ceil( (double)(RO+overlapN*(nodeN-1))/(double)nodeN );

                //size_t numOfBytesPerJob = sizeof(T)*( E1*E2*srcCHA*dstCHA*jobN + 2*E1*E2*srcCHA*jobN );

                //while ( numOfBytesPerJob > 2.2*1024*1024*1024-64.0*1024*1024 )
                //{
                //    nodeN *= 2;
                //    jobN = std::ceil( (double)(RO+overlapN*(nodeN-1))/(double)nodeN );
                //    numOfBytesPerJob = sizeof(T)*( E1*E2*srcCHA*dstCHA*jobN + 2*E1*E2*srcCHA*jobN );
                //}

                //GDEBUG_STREAM("SPIRIT - 3DT - every job will have " << numOfBytesPerJob/1024.0/1024 << " MBytes ... ");

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

                GDEBUG_STREAM("SPIRIT - 3DT - total job : " << jobList.size() << " - job N : " << jobN << " - cloud size : " << cloudSize);

                unsigned int numOfJobRunOnCloud = (unsigned int)(jobList.size() - jobList.size()/(cloudSize+1));
                if ( !runJobsOnLocalNode ) numOfJobRunOnCloud = (unsigned int)jobList.size();

                typedef Gadgetron::GadgetCloudController< gtPlusReconJob2DT<T> > GTCloudControllerType;
                GTCloudControllerType controller;

                if (controller.open () == -1)
                {
                    GERROR_STREAM("Cloud controller cannot open the cloud ...");
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
                        GDEBUG_STREAM("--> job " << j << " runs on node " << node_ids[j] << " ... ");
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
                        GERROR_STREAM("Cloud controller creates connectors failed ...");
                        controller.handle_close (ACE_INVALID_HANDLE, 0);
                        runJobsOnCloud = false;
                    }
                    else if ( controller.connectToCloud(workOrder3DT->gt_cloud_) != 0 )
                    {
                        GERROR_STREAM("Cloud controller cannot connect to the cloud ...");
                        controller.handle_close (ACE_INVALID_HANDLE, 0);
                        runJobsOnCloud = false;
                    }
                    else
                    {
                        if ( controller.runJobsOnCloud(jobListCloud, completedJobListCloud, node_ids) != 0 )
                        {
                            GERROR_STREAM("Cloud controller runs jobs on the cloud failed ...");
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
                                GDEBUG_STREAM("SPIRIT - 3DT - job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                                if ( performTiming_ ) { gt_timer3_.start("SPIRIT 3DT ... "); }
                                GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(jobList[j]));
                                if ( performTiming_ ) { gt_timer3_.stop(); }

                                std::ostringstream ostr;
                                ostr << "job_fullkspace" << "_" << j;
                                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(jobList[j].res, debugFolder_+ostr.str()); }
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
                                    GDEBUG_STREAM("SPIRIT - 3DT - uncompleted cloud job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                                    if ( performTiming_ ) { gt_timer3_.start("SPIRIT 3DT ... "); }
                                    GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(jobList[j]));
                                    if ( performTiming_ ) { gt_timer3_.stop(); }

                                    std::ostringstream ostr;
                                    ostr << "job_fullkspace" << "_" << j;
                                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(jobList[j].res, debugFolder_+ostr.str()); }
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

                GDEBUG_STREAM("SPIRIT - 3DT - total job : " << jobList.size());

                size_t j;
                for ( j=0; j<jobList.size(); j++ )
                {
                    GDEBUG_STREAM("SPIRIT - 3DT - job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                    if ( performTiming_ ) { gt_timer3_.start("SPIRIT 3DT ... "); }
                    GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(jobList[j]));
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    std::ostringstream ostr;
                    ostr << "job_fullkspace" << "_" << j;
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(jobList[j].res, debugFolder_+ostr.str()); }
                }
            }

            // combine the job
            workOrder3DTJobSplit.fullkspace_.create(E1, E2, dstCHA, RO, N);
            GADGET_CHECK_RETURN_FALSE(this->combineReconJob(&workOrder3DTJobSplit, jobList, RO, N));

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DTJobSplit.fullkspace_, debugFolder_+"job_combined_fullkspace"); }

            // clear the memory
            jobList.clear();

            // permute the unwrapped kspace
            if ( performTiming_ ) { gt_timer3_.start("permtue RO to 1st dimension ... "); }

            {
                std::vector<size_t> dim_order(4);
                dim_order[0] = 3;
                dim_order[1] = 0;
                dim_order[2] = 1;
                dim_order[3] = 2;

                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&workOrder3DTJobSplit.fullkspace_, &kspaceIfftRO, &dim_order));
            }

            if ( performTiming_ ) { gt_timer3_.stop(); }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIfftRO, debugFolder_+"res_fullkspace_ROinIm"); }

            // perform fft along the first dimension
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(kspaceIfftRO, workOrder3DT->fullkspace_);
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->fullkspace_, debugFolder_+"res_3DSpirit"); }
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

                if ( performTiming_ ) { gt_timer3_.start("spirit 3D unwrapping ... "); }
                GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(workOrder3DT, aliasedKSpace, kIm, unwarppedKSpace, n));
                if ( performTiming_ ) { gt_timer3_.stop(); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(unwarppedKSpace, debugFolder_+"unwarppedKSpace"); }
            }
        }

        hoNDArray<T> complexImMultiChannel(RO, E1, E2, dstCHA, N);

        if ( (workOrder3DT->coilMap_->get_size(0)==RO) 
            && (workOrder3DT->coilMap_->get_size(1)==E1) 
            && (workOrder3DT->coilMap_->get_size(2)==E2) 
            && (workOrder3DT->coilMap_->get_size(3)==dstCHA) )
        {
            hoNDArray<T> complexImMultiChannel(RO, E1, E2, dstCHA, N);
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->fullkspace_, complexImMultiChannel);
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(complexImMultiChannel, debugFolder_+"unwarppedComplexIm"); }

            if ( performTiming_ ) { gt_timer3_.start("spirit 3D coil combination ... "); }
            // gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(complexImMultiChannel, *workOrder3DT->coilMap_, workOrder3DT->complexIm_);
            Gadgetron::coil_combine(complexImMultiChannel, *workOrder3DT->coilMap_, 3, workOrder3DT->complexIm_);
            if ( performTiming_ ) { gt_timer3_.stop(); }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"combined"); }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTSPIRIT<T>::performUnwrapping(gtPlusReconWorkOrder3DT<T>* workOrder3DT, const hoNDArray<T>& data) ... ");
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
            GDEBUG_STREAM("gtPlusReconWorker3DTSPIRIT, kerE1!=E1 || kerE2!=E2, kernel needs to be converted along E1 and E2 ... ");

            if ( performTiming_ ) { gt_timer3_.start("kernel conversion along E1 and E2 ... "); }

            kerImE1E2RO.create(E1, E2, srcCHA, dstCHA, RO, kerN);
            Gadgetron::clear(kerImE1E2RO);

            GADGET_CHECK_RETURN_FALSE(spirit_.imageDomainKernelE1E2RO(ker, (int)E1, (int)E2, kerImE1E2RO));
            kerIm = &kerImE1E2RO;

            if ( performTiming_ ) { gt_timer3_.stop(); }
        }

        res.create(kspace.get_dimensions());

        long long NUM = (long long)RO;

        #ifdef USE_OMP
            int numThreads = NUM;

            int numOpenMPProcs = omp_get_num_procs();
            GDEBUG_STREAM("gtPlusReconWorker3DTSPIRIT, numOpenMPProcs : " << numOpenMPProcs);

            int maxOpenMPThreads = omp_get_max_threads();
            GDEBUG_STREAM("gtPlusReconWorker3DTSPIRIT, maxOpenMPThreads : " << maxOpenMPThreads);

            int allowOpenMPNested = omp_get_nested();

            if ( NUM < numOpenMPProcs-2 )
            {
                omp_set_nested(1);
                allowOpenMPNested = 1;
            }
            else
            {
                omp_set_nested(0);
                allowOpenMPNested = 0;
            }

            GDEBUG_STREAM("gtPlusReconWorker3DTSPIRIT, allowOpenMPNested : " << allowOpenMPNested);
            GDEBUG_STREAM("gtPlusReconWorker3DTSPIRIT, numThreads : " << numThreads);

            if ( numThreads > numOpenMPProcs ) numThreads = numOpenMPProcs;
            GDEBUG_STREAM("gtPlusReconWorker3DTSPIRIT, numThreads : " << numThreads);

        #endif

        long long t;

        hoNDArray<T> ker_Shifted(kerIm);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(*kerIm, ker_Shifted);
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ker_Shifted, debugFolder_+"ker_Shifted"); }

        hoNDArray<T> kspace_Shifted(kspace);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(kspace, kspace_Shifted);
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace_Shifted, debugFolder_+"kspace_Shifted"); }

        #pragma omp parallel default(none) private(t) shared(RO, E1, E2, srcCHA, dstCHA, workOrder3DT, NUM, kspace_Shifted, ker_Shifted, res) if ( NUM > 1 ) num_threads( numThreads )
        {
            std::vector<size_t> dim(3, 1);
            dim[0] = E1;
            dim[1] = E2;
            dim[2] = srcCHA;

            boost::shared_ptr< hoSPIRIT2DOperator<T> > oper(new hoSPIRIT2DOperator<T>(&dim));
            hoSPIRIT2DOperator<T>& spirit = *oper;

            spirit.use_non_centered_fft_ = true;
            spirit.no_null_space_ = false;

            hoNDArray<T> x0(E1, E2, srcCHA);
            Gadgetron::clear(x0);

            hoLSQRSolver< hoNDArray<T> > cgSolver;
            cgSolver.set_tc_tolerance((value_type)workOrder3DT->spirit_iter_thres_);
            cgSolver.set_max_iterations(workOrder3DT->spirit_iter_max_);
            cgSolver.set_verbose(workOrder3DT->spirit_print_iter_);
            cgSolver.set_encoding_operator(oper);

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

                spirit.set_forward_kernel(*kerCurr, false);

                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(E1, E2, srcCHA, kspaceCurr.begin()));
                spirit.set_acquired_points(*acq);

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*kerCurr, debugFolder_+"spirit3D_ker"); }
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*acq, debugFolder_+"spirit3D_kspace"); }

                cgSolver.set_x0(acq.get());

                // compute rhs
                spirit.compute_righ_hand_side(*acq, b);

                // solve
                cgSolver.solve(&resCurr, &b);

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(resCurr, debugFolder_+"unwarppedKSpace_t"); }

                // restore the acquired points
                spirit.restore_acquired_kspace(*acq, resCurr);

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(resCurr, debugFolder_+"unwarppedKSpace_t_setAcq"); }
            }
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"res_Shifted"); }

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fftshift2D(res, kspace_Shifted);
        res = kspace_Shifted;

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"resPermuted"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTSPIRIT<T>::performUnwarppingImplROPermuted(...) ... ");
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
        hoNDArray<T> resDecoupled(E1, E2, dstCHA, RO);

        hoNDArray<T> kspaceIfftRO(RO, E1, E2, srcCHA);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(kspace, kspaceIfftRO);
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIfftRO, debugFolder_+"kspaceIfftRO"); }

        hoNDArray<T> kspaceIfftROPermuted(E1, E2, srcCHA, RO);

        if ( performTiming_ ) { gt_timer3_.start("permtue RO to 4th dimension ... "); }

        std::vector<size_t> dim_order(4);
        dim_order[0] = 1;
        dim_order[1] = 2;
        dim_order[2] = 3;
        dim_order[3] = 0;

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&kspaceIfftRO, &kspaceIfftROPermuted, &dim_order));

        if ( performTiming_ ) { gt_timer3_.stop(); }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIfftROPermuted, debugFolder_+"kspaceIfftROPermuted"); }

        T* pKspaceIfftROPermuted = kspaceIfftROPermuted.begin();

        T* pG_I = adj_forward_G_I.begin();

        long long NUM = (long long)RO;

        long long t;

        #pragma omp parallel default(none) private(t) shared(RO, E1, E2, srcCHA, dstCHA, workOrder3DT, NUM, resDecoupled, pKspaceIfftROPermuted, pG_I) if ( NUM > 6 ) num_threads( (int)((NUM<16) ? NUM : 16) )
        {
            hoNDArray<T> adjForG_I_Decoupled(E1, E2, srcCHA, dstCHA);
            T* pDecoupledG_I = adjForG_I_Decoupled.begin();

            std::vector<size_t> dim(3, 1);
            dim[0] = E1;
            dim[1] = E2;
            dim[2] = srcCHA;

            boost::shared_ptr< hoSPIRIT2DOperator<T> > oper(new hoSPIRIT2DOperator<T>(&dim));
            hoSPIRIT2DOperator<T>& spirit = *oper;
            spirit.no_null_space_ = false;

            hoNDArray<T> x0(E1, E2, srcCHA);
            Gadgetron::clear(x0);

            hoLSQRSolver< hoNDArray<T> > cgSolver;
            cgSolver.set_tc_tolerance((value_type)workOrder3DT->spirit_iter_thres_);
            cgSolver.set_max_iterations(workOrder3DT->spirit_iter_max_);
            cgSolver.set_verbose(workOrder3DT->spirit_print_iter_);
            cgSolver.set_encoding_operator(oper);

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

                spirit.set_forward_kernel(*ker, false);

                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(E1, E2, srcCHA, kspace_DeDecoupled.begin()));
                spirit.set_acquired_points(*acq);

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*ker, debugFolder_+"spirit3D_ker"); }
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*acq, debugFolder_+"spirit3D_kspace"); }

                cgSolver.set_x0(acq.get());

                // compute rhs
                spirit.compute_righ_hand_side(*acq, b);

                // solve
                cgSolver.solve(&resCurr, &b);

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(resCurr, debugFolder_+"unwarppedKSpace_t"); }

                // restore the acquired points
                spirit.restore_acquired_kspace(*acq, resCurr);

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(resCurr, debugFolder_+"unwarppedKSpace_t_setAcq"); }
            }
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(resDecoupled, debugFolder_+"resDecoupled"); }

        // permute the unwrapped kspace
        if ( performTiming_ ) { gt_timer3_.start("permtue RO to 1st dimension ... "); }

        {
            std::vector<size_t> dim_order(4);
            dim_order[0] = 3;
            dim_order[1] = 0;
            dim_order[2] = 1;
            dim_order[3] = 2;

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&resDecoupled, &kspaceIfftRO, &dim_order));
        }

        if ( performTiming_ ) { gt_timer3_.stop(); }

        // perform fft along the first dimension
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(kspaceIfftRO, res);
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"res_3DSpirit"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTSPIRIT<T>::performUnwarppingImpl(gtPlusReconWorkOrder3DT<T>* workOrder3DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res) ... ");
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
        GERROR_STREAM("Errors in gtPlusReconWorker3DTSPIRIT<T>::performUnwarppingImpl(job) ... ");
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

        // call the BaseClass
        GADGET_CHECK_RETURN_FALSE(BaseClass::performRecon(workOrder3DT));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DTSPIRIT<T>::performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT) ... ");
        return false;
    }

    return true;
}

}}
