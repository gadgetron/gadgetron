/** \file   gtPlusISMRMRDReconWorker.h
    \brief  Define the base class for the GtPlus worker for reconstruction
    \author Hui Xue
*/

#pragma once

#include "ismrmrd.h"

#include <string>
#include "util/gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorkOrder.h"
#include "gtPlusMemoryManager.h"
#include "hoNDArrayMemoryManaged.h"
#include "SerializableObject.h"
#include "gtPlusCloudScheduler.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron { namespace gtPlus {

template <typename T> 
struct gtPlusReconJob2DT : public SerializableObject
{
    gtPlusReconWorkOrder<T> workOrder2DT;
    hoNDArray<T> kspace;
    hoNDArray<T> ker;
    // hoNDArray<T> coilMap;

    hoNDArray<T> complexIm;
    hoNDArray<T> res;

    size_t job_index_startN_;
    size_t job_index_endN_;
    size_t job_index_S_;

    gtPlusReconJob2DT();
    gtPlusReconJob2DT(const gtPlusReconJob2DT& job);

    ~gtPlusReconJob2DT();

    virtual bool serialize(char*& buf, size_t& len);
    virtual bool deserialize(char* buf, size_t& len);
};

template <typename T> 
gtPlusReconJob2DT<T>::gtPlusReconJob2DT()
{

}

template <typename T> 
gtPlusReconJob2DT<T>::~gtPlusReconJob2DT()
{

}

template <typename T> 
gtPlusReconJob2DT<T>::gtPlusReconJob2DT(const gtPlusReconJob2DT& job)
{
    job.workOrder2DT.duplicate(workOrder2DT);
    workOrder2DT.coilMap_ = job.workOrder2DT.coilMap_;
    kspace = job.kspace;
    ker = job.ker;
    // coilMap = job.coilMap;
    complexIm = job.complexIm;
    res = job.res;
    job_index_startN_ = job.job_index_startN_;
    job_index_endN_ = job.job_index_endN_;
    job_index_S_ = job.job_index_S_;
}

template <typename T> 
bool gtPlusReconJob2DT<T>::serialize(char*& buf, size_t& len)
{
    char *bufKSpace(NULL), *bufKernel(NULL), *bufCoilMap(NULL), *bufComplexIm(NULL), *bufRes(NULL);
    try
    {
        if ( buf != NULL ) delete[] buf;

        // find the total len
        gtPlusReconWorkOrderPara para;
        para = this->workOrder2DT;

        // buffer for kspace, kernel and coil map
        size_t lenKSpace, lenKernel, lenCoilMap, lenComplexIm, lenRes;

        GADGET_CHECK_THROW(kspace.serialize(bufKSpace, lenKSpace));
        GADGET_CHECK_THROW(ker.serialize(bufKernel, lenKernel));

        if ( workOrder2DT.coilMap_ )
        {
            GADGET_CHECK_THROW(workOrder2DT.coilMap_->serialize(bufCoilMap, lenCoilMap));
        }
        else
        {
            hoNDArray<T> coilMapDummy;
            GADGET_CHECK_THROW(coilMapDummy.serialize(bufCoilMap, lenCoilMap));
        }
        GADGET_CHECK_THROW(complexIm.serialize(bufComplexIm, lenComplexIm));
        GADGET_CHECK_THROW(res.serialize(bufRes, lenRes));

        // total length
        len = sizeof(gtPlusReconWorkOrderPara) + sizeof(size_t)*3 + lenKSpace + lenKernel + lenCoilMap + lenComplexIm + lenRes;

        buf = new char[len];
        GADGET_CHECK_RETURN_FALSE( buf != NULL );

        size_t offset = 0, currLen=0;

        currLen = sizeof(gtPlusReconWorkOrderPara);
        memcpy(buf+offset, &para, currLen);
        offset += currLen;

        currLen = sizeof(size_t);
        memcpy(buf+offset, &job_index_startN_, currLen);
        offset += currLen;

        currLen = sizeof(size_t);
        memcpy(buf+offset, &job_index_endN_, currLen);
        offset += currLen;

        currLen = sizeof(size_t);
        memcpy(buf+offset, &job_index_S_, currLen);
        offset += currLen;

        currLen = lenKSpace;
        memcpy(buf+offset, bufKSpace, currLen);
        offset += currLen;
        delete [] bufKSpace;

        currLen = lenKernel;
        memcpy(buf+offset, bufKernel, currLen);
        offset += currLen;
        delete [] bufKernel;

        currLen = lenCoilMap;
        memcpy(buf+offset, bufCoilMap, currLen);
        offset += currLen;
        delete [] bufCoilMap;

        currLen = lenComplexIm;
        memcpy(buf+offset, bufComplexIm, currLen);
        offset += currLen;
        delete [] bufComplexIm;

        currLen = lenRes;
        memcpy(buf+offset, bufRes, currLen);
        offset += currLen;
        delete [] bufRes;
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors happened in gtPlusReconJob2DT<T>::serialize(...) ... ");

        if ( bufKSpace != NULL ) delete [] bufKSpace;
        if ( bufKernel != NULL ) delete [] bufKernel;
        if ( bufCoilMap != NULL ) delete [] bufCoilMap;
        if ( bufComplexIm != NULL ) delete [] bufComplexIm;
        if ( bufRes != NULL ) delete [] bufRes;

        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconJob2DT<T>::deserialize(char* buf, size_t& len)
{
    try
    {
        gtPlusReconWorkOrderPara para;
        memcpy(&para, buf, sizeof(gtPlusReconWorkOrderPara));

        workOrder2DT.copyFromPara(para);

        size_t offset(sizeof(gtPlusReconWorkOrderPara)), currLen=0;

        currLen = sizeof(size_t);
        memcpy(&job_index_startN_, buf+offset, currLen);
        offset += currLen;

        currLen = sizeof(size_t);
        memcpy(&job_index_endN_, buf+offset, currLen);
        offset += currLen;

        currLen = sizeof(size_t);
        memcpy(&job_index_S_, buf+offset, currLen);
        offset += currLen;

        // kspace, kernel and coil map
        GADGET_CHECK_RETURN_FALSE(kspace.deserialize(buf+offset, currLen));
        offset += currLen;

        GADGET_CHECK_RETURN_FALSE(ker.deserialize(buf+offset, currLen));
        offset += currLen;

        hoNDArray<T> coilMapDummy;
        GADGET_CHECK_RETURN_FALSE(coilMapDummy.deserialize(buf+offset, currLen));
        offset += currLen;

        if ( coilMapDummy.get_number_of_elements() > 0 )
        {
            if ( workOrder2DT.coilMap_ )
            {
                *workOrder2DT.coilMap_ = coilMapDummy;
            }
            else
            {
                workOrder2DT.coilMap_ = boost::shared_ptr< hoNDArray<T> >( new hoNDArray<T>(coilMapDummy) );
            }
        }
        else
        {
            if ( workOrder2DT.coilMap_ ) workOrder2DT.coilMap_->clear();
        }

        GADGET_CHECK_RETURN_FALSE(complexIm.deserialize(buf+offset, currLen));
        offset += currLen;

        GADGET_CHECK_RETURN_FALSE(res.deserialize(buf+offset, currLen));
        offset += currLen;

        // total length
        len = offset;
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors happended in gtPlusReconJob2DT<T>::deserialize(...) ...");
        return false;
    }

    return true;
}

template <typename T> 
class gtPlusReconWorker
{
public:

    typedef typename realType<T>::Type value_type;

    gtPlusReconWorker() : performTiming_(false)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);
    }

    virtual ~gtPlusReconWorker() {}

    virtual bool performRecon(gtPlusReconWorkOrder<T>* workOrder) = 0;

    virtual bool performPartialFourierHandling(gtPlusReconWorkOrder<T>* /*workOrder*/) { return true; }

    virtual bool autoReconParameter(gtPlusReconWorkOrder<T>* workOrder)
    {
        if ( workOrder == NULL ) return false;
        return true;
    }

    // clock for timing
    Gadgetron::GadgetronTimer gt_timer1_;
    Gadgetron::GadgetronTimer gt_timer2_;
    Gadgetron::GadgetronTimer gt_timer3_;

    bool performTiming_;

    // exporter
    Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

    // debug folder
    std::string debugFolder_;

    // util
    gtPlusISMRMRDReconUtil<T> gtPlus_util_;

    // memory manager
    boost::shared_ptr<gtPlusMemoryManager> gtPlus_mem_manager_;

    // ----------------------------------------------------
    // recon job splitter and combiner
    // ----------------------------------------------------
    // 2DT array, [RO E1 CHA N S]
    // if splitByS is true, split jobs by each S
    // if jobN > 0, every jobN 2D kspaces are assigned into one job
    // if splitByS=false and jobN<=0, the jobMegaBytes is used to define the maximal size of every job 
    // overlapN: the overlap along N dimension
    virtual bool splitReconJob(gtPlusReconWorkOrder<T>* workOrder2DT, hoNDArray<T>& kspace, hoNDArray<T>& ker, 
                        bool splitByS, size_t jobN, size_t jobMegaBytes, size_t overlapN, 
                        std::vector<gtPlusReconJob2DT<T> >& jobList);

    virtual bool combineReconJob(gtPlusReconWorkOrder<T>* workOrder2DT, std::vector<gtPlusReconJob2DT<T> >& jobList, size_t N, size_t S);

    virtual bool createAReconJob(gtPlusReconWorkOrder<T>* workOrder2DT, hoNDArray<T>& kspace, hoNDArray<T>& ker, 
                            size_t startN, size_t endN, size_t indS, gtPlusReconJob2DT<T>& job);

    // from the node computing power indexes, get the effective node number for job splitting
    virtual bool computeEffectiveNodeNumberBasedOnComputingPowerIndex(gtPlusReconWorkOrder<T>* workOrder, size_t& numOfEffectiveNodes);

    // estimate the job size, given the maximal memory usage for every job
    virtual bool estimateJobSize(gtPlusReconWorkOrder<T>* workOrder, size_t maxNumOfBytesPerJob, size_t overlapBetweenJobs, size_t numOfNodes, size_t& jobSize) = 0;

    // given the number of nodes in a cloud and corresponding computing power indexes, spread the jobs on the nodes
    virtual bool scheduleJobForNodes(gtPlusReconWorkOrder<T>* workOrder2DT, size_t numOfJobs, std::vector<int>& nodeIdForJob);
};

template <typename T> 
bool gtPlusReconWorker<T>::createAReconJob(gtPlusReconWorkOrder<T>* workOrder2DT, hoNDArray<T>& kspace, hoNDArray<T>& ker, 
        size_t startN, size_t endN, size_t indS, gtPlusReconJob2DT<T>& job)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t N = kspace.get_size(3);
        size_t S = kspace.get_size(4);

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);
        size_t srcCHA = ker.get_size(2);
        size_t dstCHA = ker.get_size(3);
        size_t refN = ker.get_size(4);

        size_t jobN = endN-startN+1;

        job.kspace.create(RO, E1, srcCHA, jobN, 1);
        memcpy(job.kspace.begin(), kspace.begin()+indS*RO*E1*srcCHA*N+startN*RO*E1*srcCHA, job.kspace.get_number_of_bytes());

        if ( refN < N )
        {
            job.ker.create(kerRO, kerE1, srcCHA, dstCHA, refN, 1);
            memcpy(job.ker.begin(), ker.begin()+indS*kerRO*kerE1*srcCHA*dstCHA*refN, job.ker.get_number_of_bytes());
        }
        else
        {
            job.ker.create(kerRO, kerE1, srcCHA, dstCHA, jobN, 1, ker.begin()+indS*kerRO*kerE1*srcCHA*dstCHA*refN+startN*kerRO*kerE1*srcCHA*dstCHA);
        }

        if ( workOrder2DT->coilMap_->get_number_of_elements() > 0 )
        {
            if ( refN < N )
            {
                job.workOrder2DT.coilMap_ = boost::shared_ptr<hoNDArray<T> >(new hoNDArray<T>(RO, E1, dstCHA, workOrder2DT->coilMap_->begin()+indS*RO*E1*dstCHA*refN));
            }
            else
            {
                job.workOrder2DT.coilMap_ = boost::shared_ptr<hoNDArray<T> >(new hoNDArray<T>(RO, E1, dstCHA, jobN, workOrder2DT->coilMap_->begin()+indS*RO*E1*dstCHA*refN+startN*RO*E1*dstCHA));
            }
        }

        job.job_index_startN_ = startN;
        job.job_index_endN_ = endN;
        job.job_index_S_ = indS;
        workOrder2DT->duplicate(job.workOrder2DT);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker<T>::createAReconJob(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker<T>::splitReconJob(gtPlusReconWorkOrder<T>* workOrder2DT, hoNDArray<T>& kspace, hoNDArray<T>& ker, 
        bool splitByS, size_t jobN, size_t jobMegaBytes, size_t overlapN, 
        std::vector<gtPlusReconJob2DT<T> >& jobList)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t N = kspace.get_size(3);
        size_t S = kspace.get_size(4);

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);
        size_t srcCHA = ker.get_size(2);
        size_t dstCHA = ker.get_size(3);
        size_t refN = ker.get_size(4);

        size_t s;
        int startN, endN;

        if ( splitByS )
        {
            jobList.resize(S);
            startN = 0;
            endN = (int)N-1;
            for ( s=0; s<S; s++ )
            {
                GADGET_CHECK_RETURN_FALSE(createAReconJob(workOrder2DT, kspace, ker, startN, endN, s, jobList[s]));
            }

            return true;
        }

        if ( jobN > 0 )
        {
            if ( jobN < 2*overlapN ) jobN = 2*overlapN;
        }
        else if ( jobMegaBytes > 0 )
        {
            jobN = jobMegaBytes/(kerRO*kerE1*srcCHA*dstCHA*sizeof(T)/1024/1024);
            if ( jobN < 2*overlapN ) jobN = 2*overlapN;
        }

        jobList.clear();

        // find number of jobs
        size_t numPerN=0;
        startN = 0;
        while ( startN < N )
        {
            endN = (int)(startN+jobN+overlapN-1);
            numPerN++;

            if ( endN >= N )
            {
                endN = (int)N-1;
                break;
            }

            startN = endN-overlapN+1;
        }

        jobList.resize(S*numPerN);

        for ( s=0; s<S; s++ )
        {

            size_t num=0;
            startN = 0;
            while ( startN < N )
            {
                endN = (int)(startN+jobN+overlapN-1);
                num++;

                if ( endN >= N )
                {
                    endN = (int)N-1;

                    if ( endN-startN+1 < jobN )
                    {
                        startN = endN-(int)jobN+1;
                        if ( startN < 0 ) startN = 0;
                    }

                    GADGET_CHECK_RETURN_FALSE(createAReconJob(workOrder2DT, kspace, ker, startN, endN, s, jobList[s*numPerN+num-1]));
                    break;
                }

                GADGET_CHECK_RETURN_FALSE(createAReconJob(workOrder2DT, kspace, ker, startN, endN, s, jobList[s*numPerN+num-1]));

                startN = endN-(int)overlapN+1;
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker<T>::splitReconJob(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker<T>::
combineReconJob(gtPlusReconWorkOrder<T>* workOrder2DT, std::vector<gtPlusReconJob2DT<T> >& jobList, size_t N, size_t S)
{
    try
    {
        size_t RO = jobList[0].kspace.get_size(0);
        size_t E1 = jobList[0].kspace.get_size(1);

        size_t srcCHA = jobList[0].ker.get_size(2);
        size_t dstCHA = jobList[0].ker.get_size(3);
        size_t refN = jobList[0].ker.get_size(4);

        workOrder2DT->complexIm_.create(RO, E1, N, S);
        Gadgetron::clear(workOrder2DT->complexIm_);

        workOrder2DT->fullkspace_.create(RO, E1, dstCHA, N, S);
        Gadgetron::clear(workOrder2DT->fullkspace_);

        size_t ii, n, s;

        size_t numOfJobs = jobList.size();

        ho2DArray<T> fillingTimes(N, S);
        Gadgetron::clear(fillingTimes);

        for ( ii=0; ii<numOfJobs; ii++ )
        {
            size_t startN = jobList[ii].job_index_startN_;
            size_t endN = jobList[ii].job_index_endN_;
            size_t indS = jobList[ii].job_index_S_;

            if ( jobList[ii].complexIm.get_number_of_elements() > 0 )
            {
                hoNDArray<T> complexIm(RO, E1, endN-startN+1, workOrder2DT->complexIm_.begin()+indS*RO*E1*N+startN*RO*E1);
                Gadgetron::add(jobList[ii].complexIm, complexIm, complexIm);
            }

            if ( jobList[ii].res.get_number_of_elements() > 0 )
            {
                hoNDArray<T> fullkspace(RO, E1, dstCHA, endN-startN+1, workOrder2DT->fullkspace_.begin()+indS*RO*E1*dstCHA*N+startN*RO*E1*dstCHA);
                Gadgetron::add(jobList[ii].res, fullkspace, fullkspace);
            }

            for ( n=startN; n<=endN; n++ )
            {
                fillingTimes(n, indS) = fillingTimes(n, indS) + T(1.0);
            }
        }

        for ( s=0; s<S; s++ )
        {
            for ( n=0; n<N; n++ )
            {
                if ( fillingTimes(n, s).real() > 1 )
                {
                    hoNDArray<T> complexIm(RO, E1, workOrder2DT->complexIm_.begin()+s*RO*E1*N+n*RO*E1);
                    Gadgetron::scal( (value_type)(1.0)/fillingTimes(n, s).real(), complexIm);

                    hoNDArray<T> fullkspace(RO, E1, dstCHA, workOrder2DT->fullkspace_.begin()+s*RO*E1*dstCHA*N+n*RO*E1*dstCHA);
                    Gadgetron::scal( (value_type)(1.0)/fillingTimes(n, s).real(), fullkspace);
                }
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker<T>::combineReconJob(gtPlusReconWorkOrder<T>* workOrder2DT, std::vector<gtPlusReconJob2DT<T> >& jobList, size_t N, size_t S) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker<T>::
computeEffectiveNodeNumberBasedOnComputingPowerIndex(gtPlusReconWorkOrder<T>* workOrder, size_t& numOfEffectiveNodes)
{
    try
    {
        size_t numOfNodes = workOrder->gt_cloud_.size();
        numOfEffectiveNodes = 0;

        if ( numOfNodes == 0 )
        {
            GADGET_WARN_MSG("numOfNodes == 0");
            return true;
        }

        double minPowerIndex = workOrder->gt_cloud_[0].get<3>();
        double totalPowerIndex = minPowerIndex;

        size_t ii;
        for ( ii=1; ii<numOfNodes; ii++ )
        {
            totalPowerIndex += workOrder->gt_cloud_[ii].get<3>();
            if ( workOrder->gt_cloud_[ii].get<3>() < minPowerIndex ) minPowerIndex = workOrder->gt_cloud_[ii].get<3>();
        }

        numOfEffectiveNodes = (size_t)(std::floor(totalPowerIndex/minPowerIndex));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker<T>::computeEffectiveNodeNumberBasedOnComputingPowerIndex(gtPlusReconWorkOrder<T>* workOrder, unsigned int& numOfEffectiveNodes) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker<T>::
scheduleJobForNodes(gtPlusReconWorkOrder<T>* workOrder, size_t numOfJobs, std::vector<int>& nodeIdForJob)
{
    try
    {
        size_t numOfNodes = workOrder->gt_cloud_.size();

        gtPlusCloudScheduler scheduler;
        scheduler.setNumOfJobs(numOfJobs);

        std::vector<double> powerIndexes(numOfNodes);
        for ( size_t ii=0; ii<numOfNodes; ii++ )
        {
            powerIndexes[ii] = workOrder->gt_cloud_[ii].get<3>();
        }

        scheduler.setUpNodes(powerIndexes);

        GADGET_CHECK_RETURN_FALSE(scheduler.schedulerJobs(nodeIdForJob));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorker<T>::scheduleJobForNodes(gtPlusReconWorkOrder<T>* workOrder2DT, size_t numOfJobs, std::vector<int>& nodeIdForJob) ... ");
        return false;
    }

    return true;
}

}}
