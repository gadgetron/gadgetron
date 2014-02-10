/** \file   gtPlusISMRMRDReconWorker2DTSPIRIT.h
    \brief  Implement the 2DT linear SPIRIT reconstruction
    \author Hui Xue
*/

#pragma once

#include "ismrmrd.h"
#include "GadgetronTimer.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorker2DT.h"
#include "gtPlusSPIRIT.h"
#include "gtPlusSPIRIT2DTOperator.h"
#include "gtPlusLSQRSolver.h"

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

        workOrder2DT->kernel_->create(kRO, kE1, srcCHA, dstCHA, 1, 1, refN, S);
        workOrder2DT->kernelIm_->create(RO, E1, srcCHA, dstCHA, refN, S);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker2DTSPIRIT<T>::performCalibPrep(...) ... ");
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

        ho3DArray<T> acsSrc(refRO, refE1, srcCHA, const_cast<T*>(ref_src.begin()+n*refRO*refE1*srcCHA+usedS*refRO*refE1*srcCHA*refN));
        ho3DArray<T> acsDst(refRO, refE1, dstCHA, const_cast<T*>(ref_dst.begin()+n*refRO*refE1*dstCHA+usedS*refRO*refE1*dstCHA*refN));

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, acsSrc, "acsSrc");
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, acsDst, "acsDst");

        ho6DArray<T> ker(kRO, kE1, srcCHA, dstCHA, 1, 1, 
                            workOrder2DT->kernel_->begin()
                            +n*kRO*kE1*srcCHA*dstCHA
                            +usedS*kRO*kE1*srcCHA*dstCHA*refN);

        gtPlusSPIRIT2DOperator<T> spirit;
        spirit.setMemoryManager(gtPlus_mem_manager_);

        spirit.calib_use_gpu_ = workOrder2DT->spirit_use_gpu_;

        spirit.calib(acsSrc, acsDst, workOrder2DT->spirit_reg_lamda_, kRO, kE1, 1, 1, ker);

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, ker, "ker");

        bool minusI = true;

        hoNDArray<T> kIm(RO, E1, srcCHA, dstCHA, workOrder2DT->kernelIm_->begin()+n*RO*E1*srcCHA*dstCHA+usedS*RO*E1*srcCHA*dstCHA*refN);
        GADGET_CHECK_RETURN_FALSE(spirit.imageDomainKernel(ker, kRO, kE1, 1, 1, RO, E1, kIm, minusI));

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kIm, "kIm");
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker2DTSPIRIT<T>::performCalibImpl(...) ... ");
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

        int n;

        #ifdef USE_OMP
            int numThreads = (N<64) ? N : 64;

            int numOpenMPProcs = omp_get_num_procs();
            GADGET_MSG("gtPlusReconWorker2DTSPIRIT, numOpenMPProcs : " << numOpenMPProcs);

            int maxOpenMPThreads = omp_get_max_threads();
            GADGET_MSG("gtPlusReconWorker2DTSPIRIT, maxOpenMPThreads : " << maxOpenMPThreads);

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

            GADGET_MSG("gtPlusReconWorker2DTSPIRIT, allowOpenMPNested : " << allowOpenMPNested);
            GADGET_MSG("gtPlusReconWorker2DTSPIRIT, numThreads : " << numThreads);
        #endif

        GADGET_MSG("gtPlusReconWorker2DTSPIRIT, processing starts ... ");

        hoNDArray<T> ker_Shifted(adj_forward_G_I);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(adj_forward_G_I, ker_Shifted);
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, ker_Shifted, "ker_Shifted");

        hoNDArray<T> kspace_Shifted(kspace);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifftshift2D(kspace, kspace_Shifted);
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kspace_Shifted, "kspace_Shifted");

        #ifdef GCC_OLD_FLAG
            #pragma omp parallel default(none) private(n) shared(RO, E1, srcCHA, dstCHA, kspace_Shifted, ker_Shifted, workOrder2DT, refN, N) num_threads(numThreads)
        #else
            #pragma omp parallel default(none) private(n) shared(RO, E1, srcCHA, dstCHA, kspace_Shifted, ker_Shifted, workOrder2DT, res, refN, N) num_threads(numThreads)
        #endif
        {
            gtPlusSPIRIT2DOperator<T> spirit;
            // spirit.setMemoryManager(gtPlus_mem_manager_);
            spirit.use_symmetric_spirit_ = false;
            spirit.use_non_centered_fft_ = true;

            if ( refN == 1 )
            {
                boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(RO, E1, srcCHA, dstCHA, ker_Shifted.begin()));
                spirit.setForwardKernel(ker, false);
            }

            gtPlusLSQRSolver<hoNDArray<T>, hoNDArray<T>, gtPlusSPIRIT2DOperator<T> > cgSolver;

            cgSolver.iterMax_ = workOrder2DT->spirit_iter_max_;
            cgSolver.thres_ = workOrder2DT->spirit_iter_thres_;
            cgSolver.printIter_ = workOrder2DT->spirit_print_iter_;

            cgSolver.set(spirit);

            hoNDArray<T> b(RO, E1, srcCHA);

            #pragma omp for
            for ( n=0; n<(int)N; n++ )
            {
                hoNDArray<T> unwarppedKSpace(RO, E1, dstCHA, res.begin()+n*RO*E1*dstCHA);

                int kernelN = n;
                if ( kernelN >= refN ) kernelN = refN-1;

                boost::shared_ptr<hoNDArray<T> > acq(new hoNDArray<T>(RO, E1, srcCHA, kspace_Shifted.begin()+n*RO*E1*srcCHA));
                spirit.setAcquiredPoints(acq);

                cgSolver.x0_ = acq.get();

                if ( refN > 1 )
                {
                    boost::shared_ptr<hoNDArray<T> > ker(new hoNDArray<T>(RO, E1, srcCHA, dstCHA, ker_Shifted.begin()+kernelN*RO*E1*srcCHA*dstCHA));
                    spirit.setForwardKernel(ker, false);

                    // compute rhs
                    spirit.computeRighHandSide(*acq, b);

                    // solve
                    cgSolver.solve(b, unwarppedKSpace);
                }
                else
                {
                    // compute rhs
                    spirit.computeRighHandSide(*acq, b);

                    // solve
                    cgSolver.solve(b, unwarppedKSpace);
                }

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, unwarppedKSpace, "unwarppedKSpace_n");

                // restore the acquired points
                spirit.restoreAcquiredKSpace(*acq, unwarppedKSpace);

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, unwarppedKSpace, "unwarppedKSpace_n_setAcq");
            }
        }

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res, "res_Shifted");

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fftshift2D(res, kspace_Shifted);
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kspace_Shifted, "res");
        res = kspace_Shifted;
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker2DTSPIRIT<T>::performUnwarppingImpl(...) ... ");
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
            scaleFactor /= (numOfNForScaling*std::sqrt(double(srcCHA)));
        }
        else
        {
            Gadgetron::norm2(data_dst, scaleFactor);
            scaleFactor /= (N*std::sqrt(double(srcCHA)));
        }

        workOrder2DT->spirit_ncg_scale_factor_ = scaleFactor;

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
                GADGET_MSG("SPIRIT - 2DT - size of largest job : " << jobN);
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

                //GADGET_MSG("SPIRIT - 2DT - cloudSize is " << cloudSize << " - N is " << N << " ... ");
                //unsigned int nodeN = cloudSize;
                //if ( runJobsOnLocalNode ) nodeN++;
                //GADGET_MSG("SPIRIT - 2DT - runJobsOnLocalNode is " << runJobsOnLocalNode << " - nodeN is " << nodeN << " - overlapN is " << overlapN << " ... ");

                //// adjust jobN according to cloud size
                //jobN = std::ceil( (double)(N+overlapN*(nodeN-1))/(double)nodeN );

                //size_t numOfBytesPerJob = sizeof(T)*( RO*E1*srcCHA*dstCHA*jobN + 2*RO*E1*srcCHA*jobN );

                //while ( numOfBytesPerJob > 2.0*1024*1024*1024-64.0*1024*1024 )
                //{
                //    nodeN *= 2;
                //    jobN = std::ceil( (double)N/nodeN + (double)(overlapN*(nodeN-1))/nodeN );
                //    numOfBytesPerJob = sizeof(T)*( RO*E1*srcCHA*dstCHA*jobN + 2*RO*E1*srcCHA*jobN );
                //}

                //GADGET_MSG("SPIRIT - 2DT - jobN is " << jobN << "; every job has " << numOfBytesPerJob/1024.0/1024 << " MBytes ... ");

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

                GADGET_MSG("SPIRIT - 2DT - total job : " << jobList.size() << " - job N : " << jobN << " - cloud size : " << cloudSize);

                unsigned int numOfJobRunOnCloud = jobList.size() - jobList.size()/(cloudSize+1);
                if ( !runJobsOnLocalNode ) numOfJobRunOnCloud = jobList.size();
                GADGET_MSG("SPIRIT - 2DT - numOfJobRunOnCloud : " << numOfJobRunOnCloud << " ... ");

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

                    GADGET_CHECK_RETURN_FALSE(this->scheduleJobForNodes(workOrder2DT, numOfJobRunOnCloud, node_ids));

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

                    if ( controller.createConnector(workOrder2DT->gt_cloud_, GADGET_MESSAGE_CLOUD_JOB, readers, GADGET_MESSAGE_CLOUD_JOB, writers) != 0 )
                    {
                        GADGET_ERROR_MSG("Cloud controller creates connectors failed ...");
                        controller.handle_close (ACE_INVALID_HANDLE, 0);
                        runJobsOnCloud = false;
                    }
                    else if ( controller.connectToCloud(workOrder2DT->gt_cloud_) != 0 )
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
                                GADGET_MSG("SPIRIT - 2DT - job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                                GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("SPIRIT 2DT ... "));
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
                                    GADGET_MSG("SPIRIT - 2DT - uncompleted cloud job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                                    GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("SPIRIT 3DT ... "));
                                    GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(jobList[j]));
                                    GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

                                    std::ostringstream ostr;
                                    ostr << "job_fullkspace" << "_" << j;
                                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, jobList[j].res, ostr.str());
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

                GADGET_MSG("SPIRIT - 2DT - total job : " << jobList.size());

                size_t j;
                for ( j=0; j<jobList.size(); j++ )
                {
                    GADGET_MSG("SPIRIT - 2DT - job : " << j << " - size :" << jobList[j].job_index_endN_-jobList[j].job_index_startN_+1);

                    GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("L1 SPIRIT NCG 2DT ... "));
                    GADGET_CHECK_RETURN_FALSE(this->performUnwarppingImpl(jobList[j]));
                    GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, jobList[j].res, "job_fullkspace");
                }

                // combine the job
                GADGET_CHECK_RETURN_FALSE(this->combineReconJob(workOrder2DT, jobList, N, S));

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder2DT->fullkspace_, "fullkspace");

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

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, unwarppedKSpace, "unwarppedKSpace");
            }
        }

        hoNDArrayMemoryManaged<T> complexImMultiChannel(RO, E1, dstCHA, N, gtPlus_mem_manager_);

        // perform coil combination
        for ( usedS=0; usedS<S; usedS++ )
        {
            hoNDArray<T> unwarppedKSpace(RO, E1, dstCHA, N, workOrder2DT->fullkspace_.begin()+usedS*RO*E1*dstCHA*N);

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(unwarppedKSpace, complexImMultiChannel);

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, complexImMultiChannel, "unwarppedComplexIm");

            hoNDArray<T> combined(RO, E1, N, workOrder2DT->complexIm_.begin()+usedS*RO*E1*N);

            if ( refN == N )
            {
                hoNDArray<T> coilMap(RO, E1, dstCHA, refN, workOrder2DT->coilMap_->begin()+usedS*RO*E1*dstCHA*refN);
                gtPlusISMRMRDReconUtilComplex<T>().coilCombine(complexImMultiChannel, coilMap, combined);
            }
            else
            {
                hoNDArray<T> coilMap(RO, E1, dstCHA, workOrder2DT->coilMap_->begin()+usedS*RO*E1*dstCHA*refN);
                gtPlusISMRMRDReconUtilComplex<T>().coilCombine(complexImMultiChannel, coilMap, combined);
            }

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, combined, "combined");
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker2DTSPIRIT<T>::performUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& data) ... ");
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
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker2DTSPIRIT<T>::performUnwarppingImpl(...) ... ");
        return false;
    }

    return true;
}

}}
