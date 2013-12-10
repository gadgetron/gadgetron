#ifndef GADGETSTREAMCONTROLLER_H
#define GADGETSTREAMCONTROLLER_H

#include "ace/Log_Msg.h"
#include "ace/Reactor.h"
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
#include "GadgetronConnector.h"
#include "GadgetImageMessageReader.h"
#include "GadgetImageMessageWriter.h"

typedef ACE_Module<ACE_MT_SYNCH> GadgetModule;

namespace Gadgetron{

class EXPORTGADGETTOOLS GadgetStreamController 
    : public ACE_Svc_Handler<ACE_SOCK_STREAM, ACE_MT_SYNCH>
{
public:
    GadgetStreamController()
        : stream_configured_(false)
        , notifier_ (0, this, ACE_Event_Handler::WRITE_MASK)
        , writer_task_(&this->peer())
    { }

    virtual ~GadgetStreamController()
    { 
        //ACE_DEBUG( (LM_INFO, ACE_TEXT("~GadgetStreamController() called\n")) );
    }

    //ACE_SOCK_Stream &peer (void) { return this->sock_; }

    int open (void);

    /*
    virtual ACE_HANDLE get_handle (void) const { 
    return this->sock_.get_handle (); 
    }
    */

    virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE);
    //virtual int handle_output (ACE_HANDLE fd = ACE_INVALID_HANDLE);
    virtual int handle_close (ACE_HANDLE handle,
        ACE_Reactor_Mask close_mask);

    virtual int output_ready(ACE_Message_Block* mb);

    virtual Gadget* find_gadget(std::string gadget_name);

private:
    ACE_Stream<ACE_MT_SYNCH> stream_;
    bool stream_configured_;
    WriterTask writer_task_;

    ACE_Reactor_Notification_Strategy notifier_;

    GadgetMessageReaderContainer readers_;

    std::vector<ACE_DLL_Handle*> dll_handles_;

    virtual int configure(std::string config_xml_string);
    virtual int configure_from_file(std::string config_xml_filename);

    virtual GadgetModule * create_gadget_module(const char* DLL, const char* gadget, const char* gadget_module_name);

    template <class T>  T* load_dll_component(const char* DLL, const char* component_name);

};

//template<typename JobType> 
//class GadgetCloudController : public ACE_Task<ACE_MT_SYNCH>
//{
//public:
//
//    typedef boost::tuple<std::string, std::string, std::string> CloudNodeType;
//    typedef std::vector<CloudNodeType> CloudType;
//
//    GadgetCloudController();
//    virtual ~GadgetCloudController();
//
//    // this GadgetCloudController runs in the passive mode
//    virtual int open(void* = 0);
//
//    virtual int close(unsigned long flags);
//
//    // create connector and register the reader and writer for every connector
//    int createConnector(const CloudType& cloud, 
//        size_t msgID_reader, std::vector<GadgetMessageReader*>& readers, 
//        size_t msgID_writer, std::vector<GadgetMessageWriter*>& writers);
//
//    // connect to the cloud host, need to call createConnector first
//    // hostnames: the host name or IP addresses for every node
//    // port_nos: port number for every node
//    // xmlfiles: the xml configuration file name sent to every node
//    int connectToCloud(const CloudType& cloud);
//
//    // send jobs to the node and wait for jobs to be returned
//    // for every job, the node id identify which nodes to send this job
//    // after sending all jobs, this call will block until all jobs are returned
//    int runJobsOnCloud(const std::vector<int>& node_ids);
//
//    // should be called after calling runJobsOnCloud
//    int waitForJobToComplete();
//
//    // wait for all jobs to come back
//    // all returned jobs will be put into the completed_job_list_
//    // this function will not return until all jobs are returned
//    virtual int svc(void);
//
//    // list to store jobs sent to nodes
//    std::vector<JobType*> job_list_;
//    // list to store completed jobs from the nodes
//    std::vector<JobType*> completed_job_list_;
//
//private:
//
//    // connector to every node
//    // one connector for a node
//    // node id starts from 0, and increase by 1
//    std::vector<GadgetronCloudConnector<JobType>* > cloud_connectors_;
//
//    size_t cloud_msg_id_reader_;
//    size_t cloud_msg_id_writer_;
//
//    // number of available nodes in the cloud
//    unsigned int number_of_nodes_;
//
//    // node status, 0/-1 : available/unavailable
//    std::vector<int> node_status_;
//
//    // job status, 0/-1 : completed/not completed
//    std::vector<int> job_status_;
//
//    // a condition variable to wake up the caller thread
//    ACE_Thread_Mutex mutex;
//    ACE_Condition_Thread_Mutex* cond_;
//
//    ACE_Reactor gt_cloud_rector_;
//};
//
//template <typename JobType> 
//GadgetCloudController<JobType>::GadgetCloudController() : cloud_msg_id_reader_(GADGET_MESSAGE_CLOUD_JOB), cloud_msg_id_writer_(GADGET_MESSAGE_CLOUD_JOB)
//{
//    cond_ = new ACE_Condition_Thread_Mutex(mutex, "GadgetCloudController");
//}
//
//template <typename JobType> 
//GadgetCloudController<JobType>::~GadgetCloudController()
//{
//}
//
//template <typename JobType> 
//int GadgetCloudController<JobType>::open(void* p)
//{
//    ACE_TRACE(( ACE_TEXT("GadgetCloudController::open") ));
//
//    this->reactor(&gt_cloud_rector_);
//
//    //if (!this->reactor())
//    //{
//    //    ACE_DEBUG((LM_INFO, ACE_TEXT("Setting reactor")));
//    //    this->reactor(ACE_Reactor::instance());
//    //}
//
//    return this->activate( THR_NEW_LWP | THR_JOINABLE, 1 );
//}
//
//template <typename JobType> 
//int GadgetCloudController<JobType>::close(unsigned long flags)
//{
//    int rval = 0;
//    if (flags == 1)
//    {
//        ACE_Message_Block *hangup = new ACE_Message_Block();
//        hangup->msg_type( ACE_Message_Block::MB_HANGUP );
//        if (this->putq(hangup) == -1) {
//            hangup->release();
//            ACE_ERROR_RETURN( (LM_ERROR,
//                    ACE_TEXT("%p\n"),
//                    ACE_TEXT("GadgetCloudController::close, putq")),
//                    -1);
//        }
//        rval = this->wait();
//    }
//    return rval;
//}
//
//template <typename JobType> 
//int GadgetCloudController<JobType>::createConnector(const CloudType& cloud, 
//    size_t msgID_reader, std::vector<GadgetMessageReader*>& readers, 
//    size_t msgID_writer, std::vector<GadgetMessageWriter*>& writers)
//{
//    number_of_nodes_ = cloud.size();
//
//    if ( readers.size() != number_of_nodes_ ) return -1;
//    if ( writers.size() != number_of_nodes_ ) return -1;
//
//    cloud_connectors_.resize(number_of_nodes_, NULL);
//    node_status_.resize(number_of_nodes_, -1);
//
//    cloud_msg_id_reader_ = msgID_reader;
//    cloud_msg_id_writer_ = msgID_writer;
//
//    unsigned int ii;
//    for( ii=0; ii<number_of_nodes_; ii++ )
//    {
//        GadgetronCloudConnector<JobType>* con;
//        ACE_NEW_RETURN (con, GadgetronCloudConnector<JobType>, -1);
//        cloud_connectors_[ii] = con;
//
//        cloud_connectors_[ii]->register_reader(cloud_msg_id_reader_, readers[ii] );
//        cloud_connectors_[ii]->register_writer(cloud_msg_id_writer_, writers[ii] );
//
//        cloud_connectors_[ii]->set_cloud_controller(this);
//    }
//
//    return 0;
//}
//
//template <typename JobType> 
//int GadgetCloudController<JobType>::
//connectToCloud(const CloudType& cloud)
//{
//    number_of_nodes_ = cloud.size();
//    if ( cloud_connectors_.size() != number_of_nodes_ ) return -1;
//
//    unsigned int ii;
//    for( ii=0; ii<number_of_nodes_; ii++ )
//    {
//        if ( cloud_connectors_[ii] == NULL ) return -1;
//
//        // if ( cloud_connectors_[ii].open(hostnames[ii], port_nos[ii])!=0 )
//        if ( cloud_connectors_[ii]->open(cloud[ii].get<0>(), cloud[ii].get<1>())!=0 )
//        {
//            ACE_DEBUG(( LM_ERROR, ACE_TEXT("(%p) Open connection to %s:%s failed ... \n"), cloud[ii].get<0>().c_str(), cloud[ii].get<1>().c_str()));
//        }
//        else
//        {
//            node_status_[ii] = 0;
//
//            // send the xml file
//            if (cloud_connectors_[ii]->send_gadgetron_configuration_file(cloud[ii].get<2>()) != 0)
//            {
//                ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send XML configuration to the Gadgetron cloud host %s:%s \n"), cloud[ii].get<0>().c_str(), cloud[ii].get<1>().c_str()));
//                return -1;
//            }
//        }
//    }
//
//    bool hasGoodNode = false;
//    for( ii=0; ii<number_of_nodes_; ii++ )
//    {
//        if ( node_status_[ii] == 0 )
//        {
//            hasGoodNode = true;
//            break;
//        }
//    }
//
//    if ( !hasGoodNode )
//    {
//        ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to find even one good node ... \n")));
//        return -1;
//    }
//
//    return 0;
//}
//
//template <typename JobType> 
//int GadgetCloudController<JobType>::
//runJobsOnCloud(const std::vector<int>& node_ids)
//{
//    ACE_DEBUG((LM_INFO, ACE_TEXT("(%t) GadgetCloudController : into runJobsOnCloud(...) ... \n")));
//
//    if ( job_list_.empty() )
//    {
//        ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController : job list is empty ... \n")));
//        return -1;
//    }
//
//    if ( completed_job_list_.empty() )
//    {
//        ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController : completed job list is empty ... \n")));
//        return -1;
//    }
//
//    if ( job_list_.size() != completed_job_list_.size() )
//    {
//        ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController : job list size does not match ... \n")));
//        return -1;
//    }
//
//    if ( job_list_.size() != node_ids.size() )
//    {
//        ACE_DEBUG((LM_ERROR, ACE_TEXT("GadgetCloudController : job list size does not match the node id size ... \n")));
//        return -1;
//    }
//
//    std::vector<int> node_ids_used(node_ids);
//
//    unsigned int numOfJobs = job_list_.size();
//    job_status_.resize(numOfJobs, -1);
//
//    unsigned int ii;
//    for( ii=0; ii<numOfJobs; ii++ )
//    {
//        int nodeID = node_ids_used[ii];
//        if ( nodeID == -1 )
//        {
//            job_status_[ii] = 0;
//            continue;
//        }
//
//        if ( nodeID > number_of_nodes_ )
//        {
//            nodeID %= number_of_nodes_;
//        }
//
//        while ( node_status_[nodeID] < 0 )
//        {
//            nodeID--;
//            if ( nodeID == 0 ) nodeID = number_of_nodes_;
//        }
//
//        if ( nodeID != node_ids_used[ii] ) node_ids_used[ii] = nodeID;
//
//        // send job to a node
//        GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
//                new GadgetContainerMessage<GadgetMessageIdentifier>();
//
//        m1->getObjectPtr()->id = cloud_msg_id_writer_;
//
//        GadgetContainerMessage<int>* m2 =
//                new GadgetContainerMessage<int>();
//
//        *(m2->getObjectPtr()) = ii;
//
//        GadgetContainerMessage<JobType>* m3 =
//                new GadgetContainerMessage<JobType>();
//
//        *(m3->getObjectPtr()) = *(job_list_[ii]);
//        m1->cont(m2);
//        m2->cont(m3);
//
//        if ( node_status_[nodeID] == 0 )
//        {
//            if (cloud_connectors_[nodeID]->putq(m1) == -1)
//            {
//                ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send job package %d on queue for node %d \n"), ii, nodeID));
//                return -1;
//            }
//            else
//            {
//                ACE_DEBUG((LM_INFO, ACE_TEXT("Send job %d to node %d ... \n"), ii, nodeID));
//            }
//        }
//    }
//
//    std::vector<bool> closeMsgSent(number_of_nodes_, false);
//    for( ii=0; ii<numOfJobs; ii++ )
//    {
//        unsigned int nodeID = node_ids_used[ii];
//
//        if ( !closeMsgSent[nodeID] )
//        {
//            closeMsgSent[nodeID] = true;
//
//            // send the close message for this node
//            GadgetContainerMessage<GadgetMessageIdentifier>* m = new GadgetContainerMessage<GadgetMessageIdentifier>();
//            m->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;
//
//            if (cloud_connectors_[nodeID]->putq(m) == -1)
//            {
//                ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send CLOSE package on queue for node %d \n"), nodeID));
//                return -1;
//            }
//        }
//    }
//
//    ACE_DEBUG((LM_INFO, ACE_TEXT("GadgetCloudController thread - all jobs sent ... \n")));
//
//    // block the caller thread
//    // cond_->wait();
//
//    // ACE_DEBUG((LM_INFO, ACE_TEXT("GadgetCloudController thread wakes up ... \n")));
//
//    return 0;
//}
//
//template <typename JobType> 
//int GadgetCloudController<JobType>::waitForJobToComplete()
//{
//    // block the caller thread
//    ACE_DEBUG((LM_INFO, ACE_TEXT("(%t) GadgetCloudController thread sleeps ... \n")));
//    // int ret = cond_->wait();
//
//    ACE_Message_Block *mb = 0;
//    ACE_Time_Value nowait (ACE_OS::gettimeofday ());
//
//    //collect a incoming package a package if we have one
//    while (this->getq (mb) != -1)
//    {
//        GadgetContainerMessage<GadgetMessageIdentifier>* mid =
//            AsContainerMessage<GadgetMessageIdentifier>(mb);
//
//        if (!mid)
//        {
//            ACE_DEBUG ((LM_ERROR, ACE_TEXT ("Invalid message on GadgetCloudController queue\n")));
//            mb->release();
//            cond_->signal();
//            return -1;
//        }
//
//        //Is this a shutdown message?
//        if (mid->getObjectPtr()->id == GADGET_MESSAGE_CLOSE)
//        {
//            cond_->signal();
//            return 0;
//        }
//
//        if (mid->getObjectPtr()->id == cloud_msg_id_reader_)
//        {
//            GadgetContainerMessage<int>* m_jobID =
//                AsContainerMessage<int>(mid->cont());
//
//            int jobID = *(m_jobID->getObjectPtr());
//
//            GadgetContainerMessage<JobType>* job =
//                AsContainerMessage<JobType>(mid->cont()->cont());
//
//            *(completed_job_list_[jobID]) = *(job->getObjectPtr());
//            job_status_[jobID] = 0;
//        }
//
//        mb->release();
//
//        // if all jobs are received, notice the caller thread
//        bool allJobProcessed = true;
//        for ( unsigned int ii=0; ii<job_status_.size(); ii++ )
//        {
//            if ( job_status_[ii] != 0 )
//            {
//                allJobProcessed = false;
//                break;
//            }
//        }
//
//        if ( allJobProcessed )
//        {
//            ACE_DEBUG ((LM_INFO, ACE_TEXT ("All jobs are completed and returned on GadgetCloudController queue\n")));
//            break;
//        }
//    }
//
//    ACE_DEBUG((LM_INFO, ACE_TEXT("(%t) GadgetCloudController thread wakes up ... \n")));
//    return 0;
//}
//
//template <typename JobType> 
//int GadgetCloudController<JobType>::svc(void)
//{
//    ACE_DEBUG((LM_INFO, ACE_TEXT("(%t) Into GadgetCloudController svc() ... \n")));
//
//    this->reactor()->owner(ACE_Thread::self ());//, &old_owner);
//
//    this->reactor()->reset_event_loop();
//
//    ACE_Time_Value initialDelay (3);
//    ACE_Time_Value interval (0,100);
//
//    //Handle the events
//    this->reactor()->run_reactor_event_loop();
//
//    //this->reactor()->owner(&old_owner);
//
//    ACE_DEBUG ((LM_INFO, ACE_TEXT ("(%P|%t) GadgetronConnector svc done...\n")));
//
//    //ACE_Message_Block *mb = 0;
//    //ACE_Time_Value nowait (ACE_OS::gettimeofday ());
//
//    ////collect a incoming package a package if we have one
//    //while (this->getq (mb) != -1)
//    //{
//    //    GadgetContainerMessage<GadgetMessageIdentifier>* mid =
//    //            AsContainerMessage<GadgetMessageIdentifier>(mb);
//
//    //    if (!mid)
//    //    {
//    //        ACE_DEBUG ((LM_ERROR, ACE_TEXT ("Invalid message on GadgetCloudController queue\n")));
//    //        mb->release();
//    //        cond_->signal();
//    //        return -1;
//    //    }
//
//    //    //Is this a shutdown message?
//    //    if (mid->getObjectPtr()->id == GADGET_MESSAGE_CLOSE)
//    //    {
//    //        cond_->signal();
//    //        return 0;
//    //    }
//
//    //    if (mid->getObjectPtr()->id == cloud_msg_id_reader_)
//    //    {
//    //        GadgetContainerMessage<int>* m_jobID =
//    //            AsContainerMessage<int>(mid->cont());
//
//    //        int jobID = *(m_jobID->getObjectPtr());
//
//    //        GadgetContainerMessage<JobType>* job =
//    //            AsContainerMessage<JobType>(mid->cont()->cont());
//
//    //        *(completed_job_list_[jobID]) = *(job->getObjectPtr());
//    //        job_status_[jobID] = 0;
//    //    }
//
//    //    mb->release();
//
//    //    // if all jobs are received, notice the caller thread
//    //    bool allJobProcessed = true;
//    //    for ( unsigned int ii=0; ii<job_status_.size(); ii++ )
//    //    {
//    //        if ( job_status_[ii] != 0 )
//    //        {
//    //            allJobProcessed = false;
//    //            break;
//    //        }
//    //    }
//
//    //    if ( allJobProcessed )
//    //    {
//    //        ACE_DEBUG ((LM_INFO, ACE_TEXT ("All jobs are completed and returned on GadgetCloudController queue\n")));
//    //        break;
//    //    }
//    //}
//
//    //// notice the caller thread
//    //ACE_DEBUG((LM_INFO, ACE_TEXT("Wake up GadgetCloudController thread ... \n")));
//
//    //cond_->signal();
//
//    return 0;
//}

}
#endif //GADGETSTREAMCONTROLLER_H
