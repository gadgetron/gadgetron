/** \file   gtPlusISMRMRDReconWorker2DTSPIRIT.h
    \brief  Implement the 2DT linear SPIRIT reconstruction
    \author Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"
#include "GadgetronTimer.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorker2DT.h"
#include "mri_core_spirit.h"
#include "hoSPIRIT2DOperator.h"
#include "hoLSQRSolver.h"

#include "GadgetCloudController.h"
#include "GadgetCloudJobMessageReadWrite.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorker2DTSPIRIT : public gtPlusReconWorker2DT<T>
{
public:

    typedef gtPlusReconWorker2DT<T> BaseClass;
    typedef typename realType<T>::Type value_type;

    gtPlusReconWorker2DTSPIRIT() : BaseClass() {}
    virtual ~gtPlusReconWorker2DTSPIRIT() {}

    virtual bool performCalibPrep(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, gtPlusReconWorkOrder2DT<T>* workOrder2DT);
    virtual bool performCalibImpl(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, gtPlusReconWorkOrder2DT<T>* workOrder2DT, size_t n, size_t usedS);

    virtual bool performUnwarppingImpl(gtPlusReconJob2DT<T>& job);
    virtual bool performUnwarppingImpl(gtPlusReconWorkOrder<T>* workOrder2DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res, size_t s);
    virtual bool performUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& data);

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

    gtPlusSPIRIT<T> spirit_;
};

template <typename T> 
bool gtPlusReconWorker2DTSPIRIT<T>::
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

        size_t kRO = workOrder2DT->spirit_kSize_RO_;
        size_t kE1 = workOrder2DT->spirit_kSize_E1_;
        size_t oRO = workOrder2DT->spirit_oSize_RO_;
        size_t oE1 = workOrder2DT->spirit_oSize_E1_;

        workOrder2DT->kernel_->create(2*kRO-1, 2*kE1-1, srcCHA, dstCHA, refN, S);
        workOrder2DT->kernelIm_->create(RO, E1, srcCHA, dstCHA, refN, S);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DTSPIRIT<T>::performCalibPrep(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DTSPIRIT<T>::autoReconParameter(gtPlusReconWorkOrder<T>* workOrder)
{
    gtPlusReconWorkOrder2DT<T>* workOrder2DT = dynamic_cast<gtPlusReconWorkOrder2DT<T>*>(workOrder);
    if ( workOrder2DT == NULL ) return false;

    double maxAcceFactor = workOrder2DT->acceFactorE1_;

    if ( maxAcceFactor>=6 )
    {
        workOrder2DT->spirit_iter_max_ = 150;
        workOrder2DT->spirit_iter_thres_ = 0.0015;
        workOrder2DT->spirit_reg_lamda_ = 0.005;
    }
    else if ( maxAcceFactor>=5 )
    {
        workOrder2DT->spirit_iter_max_ = 120;
        workOrder2DT->spirit_iter_thres_ = 0.0015;
        workOrder2DT->spirit_reg_lamda_ = 0.005;
    }
    else if ( maxAcceFactor>=4 )
    {
        workOrder2DT->spirit_iter_max_ = 100;
        workOrder2DT->spirit_iter_thres_ = 0.0015;
        workOrder2DT->spirit_reg_lamda_ = 0.005;
    }
    else if ( maxAcceFactor>=3 )
    {
        workOrder2DT->spirit_iter_max_ = 60;
        workOrder2DT->spirit_iter_thres_ = 0.0015;
        workOrder2DT->spirit_reg_lamda_ = 0.005;
    }
    else
    {
        workOrder2DT->spirit_iter_max_ = 50;
        workOrder2DT->spirit_iter_thres_ = 0.0015;
        workOrder2DT->spirit_reg_lamda_ = 0.005;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DTSPIRIT<T>::
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

        size_t kRO = workOrder2DT->spirit_kSize_RO_;
        size_t kE1 = workOrder2DT->spirit_kSize_E1_;
        size_t oRO = workOrder2DT->spirit_oSize_RO_;
        size_t oE1 = workOrder2DT->spirit_oSize_E1_;

        ho3DArray<T> acsSrc(refRO, refE1, srcCHA, const_cast<T*>(ref_src.begin()+n*refRO*refE1*srcCHA+usedS*refRO*refE1*srcCHA*refN));
        ho3DArray<T> acsDst(refRO, refE1, dstCHA, const_cast<T*>(ref_dst.begin()+n*refRO*refE1*dstCHA+usedS*refRO*refE1*dstCHA*refN));

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(acsSrc, debugFolder_+"acsSrc"); }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(acsDst, debugFolder_+"acsDst"); }

        T* pKernel = &((*workOrder2DT->kernel_)(0, 0, 0, 0, n, usedS));
        hoNDArray<T> ker(workOrder2DT->kernel_->get_size(0), workOrder2DT->kernel_->get_size(1), srcCHA, dstCHA, pKernel);

        bool minusI = true;

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::spirit2d_calib_convolution_kernel(acsSrc, acsDst, workOrder2DT->spirit_reg_lamda_, kRO, kE1, oRO, oE1, ker, minusI));
        if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(ker, debugFolder_ + "ker"); }

        hoNDArray<T> kIm(RO, E1, srcCHA, dstCHA, workOrder2DT->kernelIm_->begin() + n*RO*E1*srcCHA*dstCHA + usedS*RO*E1*srcCHA*dstCHA*refN);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE( Gadgetron::spirit2d_image_domain_kernel(ker, RO, E1, kIm) );
        if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(kIm, debugFolder_ + "kIm"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DTSPIRIT<T>::performCalibImpl(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DTSPIRIT<T>::
performUnwarppingImpl(gtPlusReconWorkOrder<T>* workOrder2DT, hoNDArray<T>& kspace, hoNDArray<T>& adj_forward_G_I, hoNDArray<T>& res, size_t s)
{
    try
    {
        size_t refN = adj_forward_G_I.get_size(4);

        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t N = kspace.get_size(3);

        size_t srcCHA = adj_forward_G_I.get_size(2);
        size_t dstCHA = adj_forward_G_I.get_size(3);

        res.create(kspace.get_dimensions());

        long long n;

        #ifdef USE_OMP
            int numThreads = (int)( (N<64) ? N : 64 );

            int numOpenMPProcs = omp_get_num_procs();
            GDEBUG_STREAM("gtPlusReconWorker2DTSPIRIT, numOpenMPProcs : " << numOpenMPProcs);

            if ( numThreads > numOpenMPProcs ) numThreads = numOpenMPProcs;

            int maxOpenMPThreads = omp_get_max_threads();

            GDEBUG_STREAM("gtPlusReconWorker2DTSPIRIT, maxOpenMPThreads : " << maxOpenMPThreads);

            int allowOpenMPNested = omp_get_nested();

            if ( N < numOpenMPProcs-2 )
            {
                omp_set_nested(1);
                allowOpenMPNested = 1;
            }
            else
            {
                omp_set_nested(0);
                allowOpenMPNested = 0;
            }

            GDEBUG_STREAM("gtPlusReconWorker2DTSPIRIT, allowOpenMPNested : " << allowOpenMPNested);
            GDEBUG_STREAM("gtPlusReconWorker2DTSPIRIT, numThreads : " << numThreads);
            GDEBUG_STREAM("gtPlusReconWorker2DTSPIRIT, maxOpenMPThreads : " << maxOpenMPThreads);
            GDEBUG_STREAM("gtPlusReconWorker2DTSPIRIT, numThreads : " << numThreads);
        #endif

        GDEBUG_STREAM("gtPlusReconWorker2DTSPIRIT, processing starts ... ");

        hoNDArray<T> ker_Shifted(adj_forward_G_I);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(adj_forward_G_I, ker_Shifted);
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ker_Shifted, debugFolder_+"ker_Shifted"); }

        hoNDArray<T> kspace_Shifted;
        kspace_Shifted = kspace;
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(kspace, kspace_Shifted);
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace_Shifted, debugFolder_+"kspace_Shifted"); }

        hoNDArray<T> kspace_initial_Shifted;
        bool hasInitial = false;
        if ( workOrder2DT->kspace_initial_.dimensions_equal(&kspace) )
        {
            kspace_initial_Shifted = workOrder2DT->kspace_initial_;
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(workOrder2DT->kspace_initial_, kspace_initial_Shifted);
            hasInitial = true;
        }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace_initial_Shifted, debugFolder_+"kspace_initial_Shifted"); }

        #pragma omp parallel default(none) private(n) shared(RO, E1, srcCHA, dstCHA, kspace, kspace_Shifted, kspace_initial_Shifted, ker_Shifted, workOrder2DT, res, refN, N, hasInitial) num_threads(numThreads)
        {
            std::vector<size_t> dim(3, 1);
            dim[0] = RO;
            dim[1] = E1;
            dim[2] = srcCHA;

            boost::shared_ptr< hoSPIRIT2DOperator<T> > oper(new hoSPIRIT2DOperator<T>(&dim));
            hoSPIRIT2DOperator<T>& spirit = *oper;
            spirit.use_non_centered_fft_ = true;
            spirit.no_null_space_ = false;

            if (refN == 1)
            {
                boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(RO, E1, srcCHA, dstCHA, ker_Shifted.begin()));
                spirit.set_forward_kernel(*ker, false);
            }

            hoLSQRSolver< hoNDArray<T> > cgSolver;
            cgSolver.set_tc_tolerance((value_type)workOrder2DT->spirit_iter_thres_ );
            cgSolver.set_max_iterations(workOrder2DT->spirit_iter_max_);
            cgSolver.set_verbose(workOrder2DT->spirit_print_iter_);
            cgSolver.set_encoding_operator(oper);

            hoNDArray<T> b(RO, E1, srcCHA);
            hoNDArray<T> unwarppedKSpace(RO, E1, dstCHA);

            #pragma omp for
            for ( n=0; n<(long long)N; n++ )
            {
                // check whether the kspace is undersampled
                bool undersampled = false;
                for ( size_t e1=0; e1<E1; e1++ )
                {
                    if ( (std::abs( kspace(RO/2, e1, srcCHA-1, n) ) == 0)
                        && (std::abs( kspace(RO/2, e1, 0, n) ) == 0) )
                    {
                        undersampled = true;
                        break;
                    }
                }

                if ( !undersampled )
                {
                    memcpy(res.begin()+n*RO*E1*dstCHA, kspace_Shifted.begin()+n*RO*E1*srcCHA, sizeof(T)*RO*E1*dstCHA);
                    continue;
                }

                long long kernelN = n;
                if ( kernelN >= (long long)refN ) kernelN = (long long)refN-1;

                boost::shared_ptr< hoNDArray<T> > acq(new hoNDArray<T>(RO, E1, srcCHA, kspace_Shifted.begin()+n*RO*E1*srcCHA));
                spirit.set_acquired_points(*acq);

                boost::shared_ptr< hoNDArray<T> > initialAcq;
                if ( hasInitial )
                {
                    initialAcq = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>(RO, E1, srcCHA, kspace_initial_Shifted.begin()+n*RO*E1*srcCHA));
                    cgSolver.set_x0(initialAcq.get());
                }
                else
                {
                    cgSolver.set_x0(acq.get());
                }

                if ( refN > 1 )
                {
                    boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(RO, E1, srcCHA, dstCHA, ker_Shifted.begin()+kernelN*RO*E1*srcCHA*dstCHA));
                    spirit.set_forward_kernel(*ker, false);

                    // compute rhs
                    spirit.compute_righ_hand_side(*acq, b);

                    // solve
                    cgSolver.solve(&unwarppedKSpace, &b);
                }
                else
                {
                    // compute rhs
                    spirit.compute_righ_hand_side(*acq, b);

                    // solve
                    cgSolver.solve(&unwarppedKSpace, &b);
                }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(unwarppedKSpace, debugFolder_+"unwarppedKSpace_n"); }

                // restore the acquired points
                spirit.restore_acquired_kspace(*acq, unwarppedKSpace);

                memcpy(res.begin()+n*RO*E1*dstCHA, unwarppedKSpace.begin(), unwarppedKSpace.get_number_of_bytes());

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(unwarppedKSpace, debugFolder_+"unwarppedKSpace_n_setAcq"); }
            }
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"res_Shifted"); }

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fftshift2D(res, kspace_Shifted);
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace_Shifted, debugFolder_+"res"); }
        res = kspace_Shifted;
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DTSPIRIT<T>::performUnwarppingImpl(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DTSPIRIT<T>::
performUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& data_dst)
{
    try
    {
        size_t RO = workOrder2DT->data_.get_size(0);
        size_t E1 = workOrder2DT->data_.get_size(1);
        size_t N = workOrder2DT->data_.get_size(3);
        size_t S = workOrder2DT->data_.get_size(4);

        size_t srcCHA = workOrder2DT->kernelIm_->get_size(2);
        size_t dstCHA = workOrder2DT->kernelIm_->get_size(3);

        size_t refN = workOrder2DT->kernelIm_->get_size(4);

        size_t usedS;

        // compute the scaling factor
        typename realType<T>::Type scaleFactor = 1.0;
        int numOfNForScaling = 100;
        if ( N > numOfNForScaling )
        {
            hoNDArray<T> kspaceForScaleFactor(RO, E1, srcCHA, numOfNForScaling, const_cast<T*>(data_dst.begin()));
            Gadgetron::norm2(kspaceForScaleFactor, scaleFactor);
            scaleFactor /= (value_type)(numOfNForScaling*std::sqrt(double(srcCHA)));
        }
        else
        {
            Gadgetron::norm2(data_dst, scaleFactor);
            scaleFactor /= (value_type)(N*std::sqrt(double(srcCHA)));
        }

        if ( workOrder2DT->spirit_ncg_scale_factor_ < 0 )
        {
            workOrder2DT->spirit_ncg_scale_factor_ = scaleFactor;
        }
        else
        {
            GDEBUG_STREAM("SPIRIT - 2DT - spirit_ncg_scale_factor_ is preset : " << workOrder2DT->spirit_ncg_scale_factor_ << " ... ");
        }

        // split the jobs
        bool splitByS = workOrder2DT->job_split_by_S_;
        size_t jobN = workOrder2DT->job_num_of_N_;
        size_t jobMegaBytes = workOrder2DT->job_max_Megabytes_;
        size_t overlapN = workOrder2DT->job_overlap_;
        size_t maxNumOfBytesPerJob = jobMegaBytes*1024*1024;

        if ( workOrder2DT->recon_algorithm_==ISMRMRD_SPIRIT )
        {
            overlapN = 0;
        }

        bool splitJobs = (splitByS==true || jobN>0);
        if ( !splitJobs )
        {
            if ( jobMegaBytes>0 )
            {
                size_t jobN = jobMegaBytes/(RO*E1*srcCHA*dstCHA*sizeof(T)/1024/1024);
                if ( jobN < N ) splitJobs = true;
                GDEBUG_STREAM("SPIRIT - 2DT - size of largest job : " << jobN);
            }
        }

        if ( !workOrder2DT->CloudComputing_ )
        {
            if ( jobN >= N ) splitJobs = false;
        }

        if ( splitJobs )
        {
            bool runJobsOnCloud = workOrder2DT->CloudComputing_;
            unsigned int cloudSize = workOrder2DT->CloudSize_;
            bool runJobsOnLocalNode = workOrder2DT->job_perform_on_control_node_;

            std::vector<gtPlusReconJob2DT<T> > jobList;

            if ( runJobsOnCloud )
            {
                unsigned int j;

                GADGET_CHECK_RETURN_FALSE(this->estimateJobSize(workOrder2DT, maxNumOfBytesPerJob, overlapN, cloudSize, jobN));

                //GDEBUG_STREAM("SPIRIT - 2DT - cloudSize is " << cloudSize << " - N is " << N << " ... ");
                //unsigned int nodeN = cloudSize;
                //if ( runJobsOnLocalNode ) nodeN++;
                //GDEBUG_STREAM("SPIRIT - 2DT - runJobsOnLocalNode is " << runJobsOnLocalNode << " - nodeN is " << nodeN << " - overlapN is " << overlapN << " ... ");

                //// adjust jobN according to cloud size
                //jobN = std::ceil( (double)(N+overlapN*(nodeN-1))/(double)nodeN );

                //size_t numOfBytesPerJob = sizeof(T)*( RO*E1*srcCHA*dstCHA*jobN + 2*RO*E1*srcCHA*jobN );

                //while ( numOfBytesPerJob > 2.0*1024*1024*1024-64.0*1024*1024 )
                //{
                //    nodeN *= 2;
                //    jobN = std::ceil( (double)N/nodeN + (double)(overlapN*(nodeN-1))/nodeN );
                //    numOfBytesPerJob = sizeof(T)*( RO*E1*srcCHA*dstCHA*jobN + 2*RO*E1*srcCHA*jobN );
                //}

                //GDEBUG_STREAM("SPIRIT - 2DT - jobN is " << jobN << "; every job has " << numOfBytesPerJob/1024.0/1024 << " MBytes ... ");

                // split the job
                GADGET_CHECK_RETURN_FALSE(this->splitReconJob(workOrder2DT, const_cast<hoNDArray<T>&>(data_dst), *(workOrder2DT->kernelIm_), splitByS, jobN, jobMegaBytes, overlapN, jobList));

                if ( runJobsOnLocalNode )
                {
                    while ( jobList.size() <= cloudSize )
                    {
                        jobN--;
                        jobList.clear();
                        GADGET_CHECK_RETURN_FALSE(this->splitReconJob(workOrder2DT, const_cast<hoNDArray<T>&>(data_dst), *(workOrder2DT->kernelIm_), splitByS, jobN, jobMegaBytes, overlapN, jobList));
                    }
                }

                std::vector<gtPlusReconJob2DT<T> > completedJobList(jobList.size());

                for ( j=0; j<jobList.size(); j++ )
                {
                    jobList[j].workOrder2DT.duplicate(completedJobList[j].workOrder2DT);
                    completedJobList[j].job_index_startN_ = jobList[j].job_index_startN_;
                    completedJobList[j].job_index_endN_ = jobList[j].job_index_endN_;
                    completedJobList[j].job_index_S_ = jobList[j].job_index_S_;
                }

                GDEBUG_STREAM("SPIRIT - 2DT - total job : " << jobList.size() << " - job N : " << jobN << " - cloud size : " << cloudSize);

                unsigned int numOfJobRunOnCloud = (unsigned int)(jobList.size() - jobList.size()/(cloudSize+1));
                if ( !runJobsOnLocalNode ) numOfJobRunOnCloud = (unsigned int)jobList.size();
                GDEBUG_STREAM("SPIRIT - 2DT - numOfJobRunOnCloud : " << numOfJobRunOnCloud << " ... ");

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

                    GADGET_CHECK_RETURN_FALSE(this->scheduleJobForNodes(workOrder2DT, numOfJobRunOnCloud, node_ids));

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

                    if ( controller.createConnector(workOrder2DT->gt_cloud_, GADGET_MESSAGE_CLOUD_JOB, readers, GADGET_MESSAGE_CLOUD_JOB, writers) != 0 )
                    {
                        GERROR_STREAM("Cloud controller creates connectors failed ...");
                        controller.handle_close (ACE_INVALID_HANDLE, 0);
                        runJobsOnCloud = false;
                    }
                    else if ( controller.connectToCloud(workOrder2DT->gt_cloud_) != 0 )
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
                                GDEBUG_STREAM("SPIRIT - 2DT - job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                                if ( performTiming_ ) { gt_timer3_.start("SPIRIT 2DT ... "); }
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
                                    GDEBUG_STREAM("SPIRIT - 2DT - uncompleted cloud job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                                    if ( performTiming_ ) { gt_timer3_.start("SPIRIT 3DT ... "); }
                                    GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(jobList[j]));
                                    if ( performTiming_ ) { gt_timer3_.stop(); }

                                    std::ostringstream ostr;
                                    ostr << "job_fullkspace" << "_" << j;
                                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(jobList[j].res, debugFolder_+ostr.str()); }
                                }
                            }

                            // combine the job
                            GADGET_CHECK_RETURN_FALSE(this->combineReconJob(workOrder2DT, jobList, N, S));

                            // clear the memory
                            jobList.clear();
                        }
                    }
                }
            }

            if ( !runJobsOnCloud )
            {
                GADGET_CHECK_RETURN_FALSE(this->splitReconJob(workOrder2DT, const_cast<hoNDArray<T>&>(data_dst), *(workOrder2DT->kernelIm_), splitByS, jobN, jobMegaBytes, overlapN, jobList));

                GDEBUG_STREAM("SPIRIT - 2DT - total job : " << jobList.size());

                size_t j;
                for ( j=0; j<jobList.size(); j++ )
                {
                    GDEBUG_STREAM("SPIRIT - 2DT - job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                    if ( performTiming_ ) { gt_timer3_.start("L1 SPIRIT NCG 2DT ... "); }
                    GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(jobList[j]));
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(jobList[j].res, debugFolder_+"job_fullkspace"); }
                }

                // combine the job
                GADGET_CHECK_RETURN_FALSE(this->combineReconJob(workOrder2DT, jobList, N, S));

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->fullkspace_, debugFolder_+"fullkspace"); }

                // clear the memory
                jobList.clear();
            }
        }
        else
        {
            workOrder2DT->complexIm_.create(RO, E1, N, S);

            // downstream coil compression is not supported here
            // kspace is always reconed
            bool recon_kspace = true;

            workOrder2DT->fullkspace_ = data_dst;

            for ( usedS=0; usedS<S; usedS++ )
            {
                hoNDArray<T> kIm(RO, E1, srcCHA, dstCHA, refN, workOrder2DT->kernelIm_->begin()+usedS*RO*E1*srcCHA*dstCHA*refN);

                hoNDArray<T> aliasedKSpace(RO, E1, srcCHA, N, const_cast<T*>(data_dst.begin())+usedS*RO*E1*srcCHA*N);

                hoNDArray<T> unwarppedKSpace(RO, E1, dstCHA, N, workOrder2DT->fullkspace_.begin()+usedS*RO*E1*dstCHA*N);

                GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(workOrder2DT, aliasedKSpace, kIm, unwarppedKSpace, usedS));

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(unwarppedKSpace, debugFolder_+"unwarppedKSpace"); }
            }
        }

        hoNDArray<T> complexImMultiChannel(RO, E1, dstCHA, N);

        // perform coil combination
        for ( usedS=0; usedS<S; usedS++ )
        {
            hoNDArray<T> unwarppedKSpace(RO, E1, dstCHA, N, workOrder2DT->fullkspace_.begin()+usedS*RO*E1*dstCHA*N);

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(unwarppedKSpace, complexImMultiChannel);

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(complexImMultiChannel, debugFolder_+"unwarppedComplexIm"); }

            hoNDArray<T> combined(RO, E1, N, workOrder2DT->complexIm_.begin()+usedS*RO*E1*N);

            if ( refN == N )
            {
                hoNDArray<T> coilMap(RO, E1, dstCHA, refN, workOrder2DT->coilMap_->begin()+usedS*RO*E1*dstCHA*refN);
                // gtPlusISMRMRDReconUtilComplex<T>().coilCombine(complexImMultiChannel, coilMap, combined);
                Gadgetron::coil_combine(complexImMultiChannel, coilMap, 2, combined);
            }
            else
            {
                hoNDArray<T> coilMap(RO, E1, dstCHA, workOrder2DT->coilMap_->begin()+usedS*RO*E1*dstCHA*refN);
                // gtPlusISMRMRDReconUtilComplex<T>().coilCombine(complexImMultiChannel, coilMap, combined);
                Gadgetron::coil_combine(complexImMultiChannel, coilMap, 2, combined);
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(combined, debugFolder_+"combined"); }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DTSPIRIT<T>::performUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& data) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DTSPIRIT<T>::
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
        GERROR_STREAM("Errors in gtPlusReconWorker2DTSPIRIT<T>::performUnwarppingImpl(...) ... ");
        return false;
    }

    return true;
}

}}
