
#pragma once

#include "ace/Log_Msg.h"
#include "ace/Synch.h"
#include "ace/Reactor.h"
#include "ace/WFMO_Reactor.h"
#include "ace/TP_Reactor.h"
#include "ace/SOCK_Stream.h"
#include "ace/Stream.h"
#include "ace/Message_Queue.h"
#include "ace/Svc_Handler.h"
#include "ace/Reactor_Notification_Strategy.h"

#include <complex>
#include <vector>
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"

#include "gadgettools_export.h"
#include "Gadgetron.h"
#include "Gadget.h"
#include "GadgetMessageInterface.h"
#include "GadgetronCloudConnector.h"
#include "GadgetImageMessageReader.h"
#include "GadgetImageMessageWriter.h"

namespace Gadgetron
{

template<typename JobType> 
class GadgetCloudJobProcessHandler
{
public:

    GadgetCloudJobProcessHandler() {}
    virtual ~GadgetCloudJobProcessHandler() {}

    virtual bool processJob(int jobID, JobType& ajob) { return true; }
};

template<typename JobType> 
class GadgetCloudController : public ACE_Task<ACE_MT_SYNCH>
{
public:

    typedef boost::tuple<std::string, std::string, std::string, unsigned int> CloudNodeType;
    typedef std::vector<CloudNodeType> CloudType;

    GadgetCloudController();
    virtual ~GadgetCloudController();

    // this GadgetCloudController runs in the passive mode
    virtual int open(void* = 0);

    virtual int close(unsigned long flags);

    // create connector and register the reader and writer for every connector
    int createConnector(const CloudType& cloud, 
        size_t msgID_reader, std::vector<GadgetMessageReader*>& readers, 
        size_t msgID_writer, std::vector<GadgetMessageWriter*>& writers);

    // connect to the cloud host, need to call createConnector first
    // hostnames: the host name or IP addresses for every node
    // port_nos: port number for every node
    // xmlfiles: the xml configuration file name sent to every node
    int connectToCloud(const CloudType& cloud);

    // send jobs to the node and wait for jobs to be returned
    // for every job, the node id identify which nodes to send this job
    // this call can be called repeated and the wait function will wait for all jobs ever sent
    int runJobsOnCloud(std::vector<JobType*>& job_list, std::vector<JobType*>& completed_job_list, const std::vector<int>& node_ids);
    // function to ease the calling
    int runJobsOnCloud(std::vector<JobType>& job_list, std::vector<JobType>& completed_job_list, const std::vector<int>& node_ids);

    // should be called after calling runJobsOnCloud
    int waitForJobToComplete();

    // send close message to all nodes
    int closeCloudNode();

    virtual int handle_close (ACE_HANDLE handle, ACE_Reactor_Mask close_mask);

    // set jobs on a node to be completed
    // if jobID===-1, all jobs for this node is set to be completed
    int setJobsTobeCompleted(unsigned int nodeID, int jobID=-1);

    // get/set the node status, 0/-1 : available/unavailable
    int get_node_status(int nodeID, int& status)
    {
        ACE_GUARD_RETURN(ACE_Thread_Mutex, guard, cloud_controller_mutex_, -1);
        if ( (nodeID>=0) && (nodeID<node_status_.size()) )
        {
            status = node_status_[nodeID];
        }
        else
        {
            status = -1;
        }
        return 0;
    }

    int set_node_status(int nodeID, int status)
    {
        ACE_GUARD_RETURN(ACE_Thread_Mutex, guard, cloud_controller_mutex_, -1);
        if ( (nodeID>=0) && (nodeID<node_status_.size()) ) node_status_[nodeID] = status;
        return 0;
    }

    // append the job list
    int appendJobList(std::vector<JobType*>& job_list, 
        std::vector<JobType*>& completed_job_list, 
        std::vector<int>& node_id_used, std::vector<int>& job_status);

    // list to store jobs sent to nodes
    std::vector<JobType*> job_list_;
    // list to store completed jobs from the nodes
    std::vector<JobType*> completed_job_list_;
    // for every job, indicate which node a job is sent to
    std::vector<int> node_id_used_;
    // job status, 0/-1 : completed/not completed
    std::vector<int> job_status_;

    // a function handler to process job after receive
    // this is a hook to give user a chance to do some processing after receiving every job
    GadgetCloudJobProcessHandler<JobType>* job_handler_;

private:

    // connector to every node
    // one connector for a node
    // node id starts from 0, and increase by 1
    std::vector<GadgetronCloudConnector<JobType>* > cloud_connectors_;

    size_t cloud_msg_id_reader_;
    size_t cloud_msg_id_writer_;

    // number of available nodes in the cloud
    unsigned int number_of_nodes_;

    // node status, 0/-1 : available/unavailable
    std::vector<int> node_status_;

    // number of job actually sent to nodes
    // if 0, then controller does not need to wait
    unsigned int number_of_jobs_sent_out_;

    // to protect the access to job_status_ and node_id_used_
    ACE_Thread_Mutex cloud_controller_mutex_;
};

template <typename JobType> 
GadgetCloudController<JobType>::GadgetCloudController() : cloud_msg_id_reader_(GADGET_MESSAGE_CLOUD_JOB), cloud_msg_id_writer_(GADGET_MESSAGE_CLOUD_JOB), job_handler_(NULL), number_of_jobs_sent_out_(0)
{

}

template <typename JobType> 
GadgetCloudController<JobType>::~GadgetCloudController()
{
    GADGET_DEBUG1("Into ~GadgetCloudController() ... \n");
    this->msg_queue()->deactivate();

    for ( unsigned int ii=0; ii<cloud_connectors_.size(); ii++ )
    {
        if ( cloud_connectors_[ii] != NULL )
        {
            cloud_connectors_[ii]->close();
            delete cloud_connectors_[ii];
            cloud_connectors_[ii] = NULL;
            GADGET_DEBUG1("~GadgetCloudController() : clean connectors done \n");
        }
    }
}

template <typename JobType> 
int GadgetCloudController<JobType>::open(void* p)
{
    GADGET_DEBUG1("GadgetCloudController::open\n");

    // set the high water mark of message queue to be 24GB
    this->msg_queue()->high_water_mark( (size_t)(24.0*1024*1024*1024) );

    return 0;
}

template <typename JobType> 
int GadgetCloudController<JobType>::close(unsigned long flags)
{
    GADGET_DEBUG1("GadgetCloudController::close\n");
    int rval = 0;
    if (flags == 1)
    {
        ACE_Message_Block *hangup = new ACE_Message_Block();
        hangup->msg_type( ACE_Message_Block::MB_HANGUP );
        if (this->putq(hangup) == -1)
        {
            hangup->release();
            ACE_ERROR_RETURN( (LM_ERROR,
                    ACE_TEXT("%p\n"),
                    ACE_TEXT("GadgetCloudController::close, putq")),
                    -1);
        }
        rval = this->wait();
    }
    return rval;
}

template <typename JobType> 
int GadgetCloudController<JobType>::createConnector(const CloudType& cloud, 
    size_t msgID_reader, std::vector<GadgetMessageReader*>& readers, 
    size_t msgID_writer, std::vector<GadgetMessageWriter*>& writers)
{
    number_of_nodes_ = (unsigned int)cloud.size();

    if ( readers.size() != number_of_nodes_ ) return -1;
    if ( writers.size() != number_of_nodes_ ) return -1;

    cloud_connectors_.resize(number_of_nodes_, NULL);
    node_status_.resize(number_of_nodes_, -1);

    cloud_msg_id_reader_ = msgID_reader;
    cloud_msg_id_writer_ = msgID_writer;

    unsigned int ii;
    for( ii=0; ii<number_of_nodes_; ii++ )
    {
        GadgetronCloudConnector<JobType>* con;
        ACE_NEW_RETURN (con, GadgetronCloudConnector<JobType>, -1);

        cloud_connectors_[ii] = con;
        cloud_connectors_[ii]->nodeID_ = ii;

        cloud_connectors_[ii]->register_reader(cloud_msg_id_reader_, readers[ii] );
        cloud_connectors_[ii]->register_writer(cloud_msg_id_writer_, writers[ii] );

        cloud_connectors_[ii]->set_cloud_controller(this);
    }

    return 0;
}

template <typename JobType> 
int GadgetCloudController<JobType>::
connectToCloud(const CloudType& cloud)
{
    number_of_nodes_ = (unsigned int)cloud.size();
    if ( cloud_connectors_.size() != number_of_nodes_ ) return -1;

    node_status_.resize(number_of_nodes_, -1);

    unsigned int ii;
    for( ii=0; ii<number_of_nodes_; ii++ )
    {
        if ( cloud_connectors_[ii] == NULL ) return -1;

        std::string host = cloud[ii].get<0>();
        std::string port = cloud[ii].get<1>();

        if ( cloud_connectors_[ii]->open(cloud[ii].get<0>(), cloud[ii].get<1>())!=0 )
        {
            cloud_connectors_[ii]->set_status(false);

            ACE_Time_Value tv( (time_t)GADGETRON_TIMEOUT_PERIOD );
            ACE_OS::sleep(tv);

            GADGET_DEBUG2("Open connection to (%s):%s failed ... \n", host.c_str(), port.c_str());
        }
        else
        {
            ACE_Time_Value tv( (time_t)0.5);
            ACE_OS::sleep(tv);

            // send the xml file
            if (cloud_connectors_[ii]->send_gadgetron_configuration_file(cloud[ii].get<2>()) != 0)
            {
                ACE_Time_Value tv( (time_t)GADGETRON_TIMEOUT_PERIOD );
                ACE_OS::sleep(tv);

                GADGET_DEBUG2("Unable to send XML configuration to the Gadgetron cloud host (%s):%s ... \n", host.c_str(), port.c_str());
            }
            else
            {
                // indicate this node can be used
                node_status_[ii] = 0;
                cloud_connectors_[ii]->set_status(true);
            }
        }

        if ( node_status_[ii] == 0 )
        {
            GADGET_DEBUG2("--> Node (%s):%s is ready ... \n", host.c_str(), port.c_str());
        }
        else
        {
            GADGET_DEBUG2("--> Node (%s):%s is NOT ready ... \n", host.c_str(), port.c_str());
        }
    }

    bool hasGoodNode = false;
    for( ii=0; ii<number_of_nodes_; ii++ )
    {
        if ( node_status_[ii] == 0 )
        {
            hasGoodNode = true;
            break;
        }
    }

    if ( !hasGoodNode )
    {
        ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to find even one good node ... \n")));
        return -1;
    }

    return 0;
}

template <typename JobType> 
int GadgetCloudController<JobType>::
runJobsOnCloud(std::vector<JobType*>& job_list, std::vector<JobType*>& completed_job_list, const std::vector<int>& node_ids)
{
    ACE_DEBUG((LM_INFO, ACE_TEXT("(%t) GadgetCloudController : into runJobsOnCloud(...) ... \n")));

    if ( job_list.empty() )
    {
        ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController : job list is empty ... \n")));
        return -1;
    }

    if ( completed_job_list.empty() )
    {
        ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController : completed job list is empty ... \n")));
        return -1;
    }

    if ( job_list.size() != completed_job_list.size() )
    {
        ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController : job list size does not match ... \n")));
        return -1;
    }

    if ( job_list.size() != node_ids.size() )
    {
        ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController : job list size does not match the node id size ... \n")));
        return -1;
    }

    std::vector<int> node_ids_used(node_ids);

    size_t numOfJobs = job_list.size();
    std::vector<int> job_status(numOfJobs, -1);

    size_t ii;
    for( ii=0; ii<numOfJobs; ii++ )
    {
        int nodeID = node_ids_used[ii];
        if ( nodeID == -1 )
        {
            job_status[ii] = 0;
            continue;
        }

        if ( nodeID >= (int)number_of_nodes_ )
        {
            nodeID %= (int)number_of_nodes_;
        }

        /*while ( node_status_[nodeID] < 0 )
        {
            nodeID--;
            if ( nodeID < 0 ) nodeID = number_of_nodes_-1;
        }

        if ( nodeID != node_ids_used[ii] ) node_ids_used[ii] = nodeID;*/

        int status = -1;
        this->get_node_status(nodeID, status);
        if ( status < 0 )
        {
            // try again
            if ( number_of_nodes_ > 1 )
            {
                nodeID += number_of_nodes_/2;
                if ( nodeID >= (int)number_of_nodes_ )
                {
                    nodeID %= (int)number_of_nodes_;
                }

                this->get_node_status(nodeID, status);
            }

            if ( status < 0 )
            {
                node_ids_used[ii] = -1; // local node to perform this job
                job_status[ii] = 0;
            }
            else
            {
                node_ids_used[ii] = nodeID;
            }
        }

        GADGET_DEBUG2("--> node for job %d is %d ... \n", ii, node_ids_used[ii]);
    }

    // append incoming jobs into the list
    size_t startJobID = job_list_.size();

    if ( this->appendJobList(job_list, completed_job_list, node_ids_used, job_status) == -1 )
    {
        ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to append job list ... \n")));
        return -1;
    }

    for( ii=0; ii<numOfJobs; ii++ )
    {
        int nodeID = node_ids_used[ii];
        if ( nodeID == -1 )
        {
            GADGET_DEBUG2("--> node for job %d is NOT ready ... \n", ii+startJobID);
            continue;
        }

        // send job to a node
        GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
                new GadgetContainerMessage<GadgetMessageIdentifier>();

        m1->getObjectPtr()->id = (ACE_UINT16)cloud_msg_id_writer_;

        GadgetContainerMessage<int>* m2 =
                new GadgetContainerMessage<int>();

        *(m2->getObjectPtr()) = (int)(ii+startJobID);

        GadgetContainerMessage<JobType>* m3 =
                new GadgetContainerMessage<JobType>();

        *(m3->getObjectPtr()) = *(job_list[ii]);
        m1->cont(m2);
        m2->cont(m3);

        if ( node_status_[nodeID] == 0 )
        {
            if (cloud_connectors_[nodeID]->putq(m1) == -1)
            {
                ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send job package %d on queue for node %d \n"), ii+startJobID, nodeID));
                m1->release();
                return -1;
            }
            else
            {
                GADGET_DEBUG2("Send job %d to node %d ... \n", ii+startJobID, nodeID);
                number_of_jobs_sent_out_++;
            }
        }
        else
        {
            m1->release();
        }
    }

    GADGET_DEBUG1("GadgetCloudController - all jobs sent ... \n");

    return 0;
}

template <typename JobType> 
int GadgetCloudController<JobType>::
runJobsOnCloud(std::vector<JobType>& job_list, std::vector<JobType>& completed_job_list, const std::vector<int>& node_ids)
{
    if ( job_list.size() != completed_job_list.size() )
    {
        ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController : job list size does not match ... \n")));
        return -1;
    }

    if ( job_list.size() != node_ids.size() )
    {
        ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController : job list size does not match the node id size ... \n")));
        return -1;
    }

    std::vector<JobType*> jobPtr(job_list.size(), NULL);
    std::vector<JobType*> completedJobPtr(completed_job_list.size(), NULL);

    unsigned int N = job_list.size();

    unsigned int ii;
    for ( ii=0; ii<N; ii++ )
    {
        jobPtr[ii] = &job_list[ii];
        completedJobPtr[ii] = &completed_job_list[ii];
    }

    return runJobsOnCloud(jobPtr, completedJobPtr, node_ids);
}

template <typename JobType> 
int GadgetCloudController<JobType>::
closeCloudNode()
{
    GADGET_DEBUG1("GadgetCloudController : into closeCloudNode(...) ... \n");

    unsigned int ii;

    std::vector<bool> closeMsgSent(number_of_nodes_, false);
    for( ii=0; ii<number_of_nodes_; ii++ )
    {
        int nodeID = ii;

        if ( !closeMsgSent[nodeID] )
        {
            closeMsgSent[nodeID] = true;

            // send the close message for this node
            GadgetContainerMessage<GadgetMessageIdentifier>* m = new GadgetContainerMessage<GadgetMessageIdentifier>();
            m->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;

            if (cloud_connectors_[nodeID]->putq(m) == -1)
            {
                ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send CLOSE package on queue for node %d \n"), nodeID));
                m->release();
                return -1;
            }
        }
    }

    GADGET_DEBUG1("GadgetCloudController - close message sent to all nodes ... \n");

    return 0;
}

template <typename JobType> 
int GadgetCloudController<JobType>::waitForJobToComplete()
{
    // block the caller thread
    GADGET_DEBUG1("GadgetCloudController waitForJobToComplete ... \n");

    ACE_Message_Block *mb = 0;
    ACE_Time_Value nowait (ACE_OS::gettimeofday ());

    //collect a incoming package a package if we have one
    while ( number_of_jobs_sent_out_>0 && (this->getq (mb) != -1) )
    {
        GadgetContainerMessage<int>* m_jobID =
            AsContainerMessage<int>(mb);

        if ( !m_jobID )
        {
            ACE_DEBUG ((LM_INFO, ACE_TEXT ("Invalid message id in the GadgetCloudController queue\n")));
            break;
        }

        int jobID = *(m_jobID->getObjectPtr());

        if ( jobID != -1 )
        {
            GadgetContainerMessage<JobType>* job =
                AsContainerMessage<JobType>(mb->cont());

            if ( !job )
            {
                ACE_DEBUG ((LM_INFO, ACE_TEXT ("Invalid message obj in the GadgetCloudController queue\n")));
                break;
            }

            *(completed_job_list_[jobID]) = *(job->getObjectPtr());
            job_status_[jobID] = 0;

            ACE_DEBUG ((LM_INFO, ACE_TEXT ("--> receive completed job : %d ... \n"), jobID));

            if ( job_handler_ != NULL )
            {
                if ( !job_handler_->processJob( jobID, *(completed_job_list_[jobID]) ) )
                {
                    ACE_DEBUG ((LM_INFO, ACE_TEXT ("job_handler_->processJob after receiving failed\n")));
                }
            }
        }
        else
        {
            ACE_DEBUG ((LM_INFO, ACE_TEXT ("--> receive jobID == -1 ... \n")));
        }

        mb->release();

        // if all jobs are received, notice the caller thread
        bool allJobProcessed = true;
        {
            ACE_GUARD_RETURN(ACE_Thread_Mutex, guard, cloud_controller_mutex_, -1);
            for ( unsigned int ii=0; ii<job_status_.size(); ii++ )
            {
                if ( job_status_[ii] != 0 )
                {
                    allJobProcessed = false;
                    break;
                }
            }
        }

        if ( allJobProcessed )
        {
            ACE_DEBUG ((LM_INFO, ACE_TEXT ("All jobs are completed and returned on GadgetCloudController queue\n")));
            break;
        }
    }

    // need to wait for all reader task to complete
    for( unsigned int ii=0; ii<number_of_nodes_; ii++ )
    {
        if ( cloud_connectors_[ii]->status() )
        {
            cloud_connectors_[ii]->wait();
        }
    }

    ACE_DEBUG((LM_INFO, ACE_TEXT("(%t) GadgetCloudController waitForJobToComplete done ... \n")));
    return 0;
}

template <typename JobType> 
int GadgetCloudController<JobType>::handle_close(ACE_HANDLE handle, ACE_Reactor_Mask close_mask)
{
    GADGET_DEBUG1("GadgetCloudController handling close...\n");
    return this->wait();
}

template<typename JobType> 
int GadgetCloudController<JobType>::setJobsTobeCompleted(unsigned int nodeID, int jobID)
{
    ACE_GUARD_RETURN(ACE_Thread_Mutex, guard, cloud_controller_mutex_, -1);
    try
    {
        if ( (nodeID>=0) && (nodeID<this->node_status_.size()) )
        {
            node_status_[nodeID] = -1;
        }

        size_t N = this->node_id_used_.size();
        size_t ii;
        for ( ii=0; ii<N; ii++ )
        {
            if ( this->node_id_used_[ii] == nodeID )
            {
                //if ( jobID>=0 && jobID<this->job_status_.size() )
                //{
                //    this->job_status_[jobID] = 0;
                //}
                //else
                //{
                //    if ( this->job_status_[ii]!= 0 ) this->job_status_[ii] = 0;
                //}

                // make sure all jobs on this node is marked as completed
                if ( this->job_status_[ii]!= 0 ) this->job_status_[ii] = 0;
            }
        }
    }
    catch(...)
    {
        ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetCloudController, setJobsTobeCompleted() failed ... \n")) );
        return -1;
    }

    return 0;
}

template<typename JobType> 
int GadgetCloudController<JobType>::appendJobList(std::vector<JobType*>& job_list, 
        std::vector<JobType*>& completed_job_list, 
        std::vector<int>& node_id_used, std::vector<int>& job_status)
{
    ACE_GUARD_RETURN(ACE_Thread_Mutex, guard, cloud_controller_mutex_, -1);
    try
    {
        size_t N = job_list.size();

        if ( completed_job_list.size() != N )
        {
            ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController appendJobList: job list size does not match ... \n")));
            return -1;
        }

        if ( node_id_used.size() != N )
        {
            ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController appendJobList: node_id_used size does not match ... \n")));
            return -1;
        }

        if ( job_status.size() != N )
        {
            ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController appendJobList: job_status size does not match ... \n")));
            return -1;
        }

        size_t ii;
        for ( ii=0; ii<N; ii++ )
        {
            job_list_.push_back(job_list[ii]);
            completed_job_list_.push_back(completed_job_list[ii]);
            node_id_used_.push_back(node_id_used[ii]);
            job_status_.push_back(job_status[ii]);
        }
    }
    catch(...)
    {
        ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetCloudController, appendJobList() failed ... \n")) );
        return -1;
    }

    return 0;
}

}
