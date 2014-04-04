/** \file   GadgetCloudJobMessageReaderWriter.h
    \brief  Implement the reader/writer for GtPlus cloud job data package
            This implementation requires the cloud job supports the serialize and deserialize functions

    \author Hui Xue
*/

#pragma once

#include "GtPlusGadgetExport.h"

#include "GadgetImageMessageReader.h"
#include "GadgetImageMessageWriter.h"
#include "gtPlusISMRMRDReconWorker.h"
#include "GtPlusRecon2DTCloudPackage.h"

namespace Gadgetron
{

    template <typename JobType> 
    class GadgetCloudJobMessageReader : public GadgetMessageReader
    {

    public:
        virtual ACE_Message_Block* read(ACE_SOCK_Stream* stream) 
        {
            GadgetContainerMessage<int>* jobID = new GadgetContainerMessage<int>();
            GadgetContainerMessage<JobType>* job = new GadgetContainerMessage<JobType>();

            jobID->cont(job);

            int id = 0;
            size_t sizeOfJob = 0;

            ssize_t recv_count = 0;

            if ((recv_count = stream->recv_n(&id, sizeof(int))) <= 0)
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetCloudJobMessageReader, failed to read job id\n")) );
                job->release();
                return 0;
            }

            *(jobID->getObjectPtr()) = id;

            if ((recv_count = stream->recv_n(&sizeOfJob, sizeof(size_t))) <= 0)
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetCloudJobMessageReader, failed to read job size\n")) );
                job->release();
                return 0;
            }

            hoNDArray<char> jobBuf;
            try
            {
                jobBuf.create(sizeOfJob);
            }
            catch(...)
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetCloudJobMessageReader, failed to allocate memory\n")) );
                job->release();
                return 0;
            }

            size_t maxBytesPerSend = 512.0*1024*1024;

            if ( sizeOfJob > maxBytesPerSend )
            {
                size_t receivedBytes = 0;
                size_t receivingBytes = maxBytesPerSend;

                while ( receivingBytes > 0 )
                {
                    if ((recv_count = stream->recv_n(jobBuf.get_data_ptr()+receivedBytes, receivingBytes)) <= 0)
                    {
                        ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetCloudJobMessageReader, failed to read data from socket\n")) );
                        job->release();
                        return 0;
                    }

                    receivedBytes += receivingBytes;
                    if ( receivedBytes >= sizeOfJob ) break;

                    if ( sizeOfJob-receivedBytes < maxBytesPerSend )
                    {
                        receivingBytes = sizeOfJob-receivedBytes;
                    }
                }
            }
            else
            {
                if ((recv_count = stream->recv_n(jobBuf.get_data_ptr(), sizeOfJob)) <= 0)
                {
                    ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetCloudJobMessageReader, failed to read data from socket\n")) );
                    job->release();
                    return 0;
                }
            }

            if ( !job->getObjectPtr()->deserialize(jobBuf.get_data_ptr(), sizeOfJob) ) return 0;

            return jobID;
        }
    };

    template <typename JobType> 
    class GadgetCloudJobMessageWriter : public GadgetMessageWriter
    {

    public:

        ACE_UINT16 msg_id_;

        GadgetCloudJobMessageWriter() : msg_id_(GADGET_MESSAGE_CLOUD_JOB) {}

        virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb) 
        {

            GadgetContainerMessage<int>* m1 = 
                dynamic_cast< GadgetContainerMessage<int>* >(mb);

            int jobID = *(m1->getObjectPtr());

            GadgetContainerMessage<JobType>* job = 
                dynamic_cast< GadgetContainerMessage<JobType>* >(mb->cont());

            if (!job )
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), GadgetCloudJobMessageWriter invalid image message objects")) );
                return -1;
            }

            ssize_t send_cnt = 0;
            GadgetMessageIdentifier id;
            id.id = msg_id_;

            if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0)
            {
                ACE_DEBUG((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send job message identifier\n")));
                return -1;
            }

            if ((send_cnt = sock->send_n (&jobID, sizeof(int))) <= 0)
            {
                ACE_DEBUG((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send job id\n")));
                return -1;
            }

            size_t sizeOfJob=0;
            char* buf = NULL;
            if ( !job->getObjectPtr()->serialize(buf, sizeOfJob) )
            {
                ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetCloudJobMessageWriter, failed to serialize the job\n")) );
                return -1;
            }

            if ((send_cnt = sock->send_n (&sizeOfJob, sizeof(size_t))) <= 0)
            {
                ACE_DEBUG ((LM_ERROR,
                    ACE_TEXT ("(%P|%t) Unable to send job size\n")));
                delete [] buf;
                return -1;
            }

            GADGET_DEBUG2("--> send job, size of job : %f MBytes ... \n", sizeOfJob/1024.0/1024);

            size_t maxBytesPerSend = 512.0*1024*1024;

            if ( sizeOfJob > maxBytesPerSend )
            {
                size_t sentBytes = 0;
                size_t sendingBytes = maxBytesPerSend;

                while ( sendingBytes > 0 )
                {
                    if ((send_cnt = sock->send_n (buf+sentBytes, sendingBytes)) <= 0)
                    {
                        ACE_DEBUG ((LM_ERROR,
                            ACE_TEXT ("(%P|%t) Unable to send job data\n")));
                        delete [] buf;
                        return -1;
                    }

                    sentBytes += sendingBytes;
                    if ( sentBytes >= sizeOfJob ) break;

                    if ( sizeOfJob-sentBytes < maxBytesPerSend )
                    {
                        sendingBytes = sizeOfJob-sentBytes;
                    }
                }
            }
            else
            {
                if ((send_cnt = sock->send_n (buf, sizeOfJob)) <= 0)
                {
                    ACE_DEBUG ((LM_ERROR,
                        ACE_TEXT ("(%P|%t) Unable to send job data\n")));
                    delete [] buf;
                    return -1;
                }
            }

            delete [] buf;

            return 0;
        }

    };

    typedef Gadgetron::gtPlus::gtPlusReconJob2DT< std::complex<float> > GtPlusReconJobTypeCPFL;

    class EXPORTGTPLUSGADGET GtPlusCloudJobMessageReaderCPFL : public GadgetCloudJobMessageReader<GtPlusReconJobTypeCPFL> {};
    class EXPORTGTPLUSGADGET GtPlusCloudJobMessageWriterCPFL : public GadgetCloudJobMessageWriter<GtPlusReconJobTypeCPFL> {};

    // gadget level cloud computing

    class EXPORTGTPLUSGADGET GtPlus2DTGadgetCloudJobMessageReaderCPFL : public GadgetCloudJobMessageReader<GtPlusRecon2DTCloudPackageCPFL> {};
    class EXPORTGTPLUSGADGET GtPlus2DTGadgetCloudJobMessageWriterCPFL : public GadgetCloudJobMessageWriter<GtPlusRecon2DTCloudPackageCPFL>
    {
    public:
        typedef GadgetCloudJobMessageWriter<GtPlusRecon2DTCloudPackageCPFL> BaseClass;

        GtPlus2DTGadgetCloudJobMessageWriterCPFL() : BaseClass()
        {
            msg_id_ = GADGET_MESSAGE_GADGETCLOUD_JOB;
        }
    };
}
