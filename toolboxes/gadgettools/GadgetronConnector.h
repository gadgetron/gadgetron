/*
 * GadgetronConnector.h
 *
 *  Created on: Nov 1, 2011
 *      Author: Michael S. Hansen
 */


#ifndef GADGETRONCONNECTOR_H_
#define GADGETRONCONNECTOR_H_

#include <ace/Svc_Handler.h>
#include <ace/Reactor.h>
#include <ace/SOCK_Stream.h>
#include <ace/Reactor_Notification_Strategy.h>


#include <string>

#include "GadgetronSlotContainer.h"
#include "GadgetMessageInterface.h"
#include "gadgettools_export.h"

class WriterTask : public ACE_Task<ACE_MT_SYNCH>
{

public:
	typedef ACE_Task<ACE_MT_SYNCH> inherited;

	WriterTask(ACE_SOCK_Stream* socket)
	: inherited()
	, socket_(socket)
	{
	}

	virtual ~WriterTask()
	{
		writers_.clear();
	}

	virtual int init(void)
	{
		ACE_TRACE(( ACE_TEXT("WriterTask::init") ));
		return 0;
	}

	virtual int open(void* = 0)
	{
		ACE_TRACE(( ACE_TEXT("WriterTask::open") ));

		return this->activate( THR_NEW_LWP | THR_JOINABLE,
				1 );
	}


	int register_writer(unsigned int slot, GadgetMessageWriter* writer) {
		return writers_.insert(slot,writer);
	}


	virtual int close(unsigned long flags)
	{
		int rval = 0;
		if (flags == 1) {
			ACE_Message_Block *hangup = new ACE_Message_Block();
			hangup->msg_type( ACE_Message_Block::MB_HANGUP );
			if (this->putq(hangup) == -1) {
				hangup->release();
				ACE_ERROR_RETURN( (LM_ERROR,
						ACE_TEXT("%p\n"),
						ACE_TEXT("WriterTask::close, putq")),
						-1);
			}
			rval = this->wait();
		}
		return rval;
	}

	virtual int svc(void)
	{
		ACE_Message_Block *mb = 0;
		ACE_Time_Value nowait (ACE_OS::gettimeofday ());


		//Send a package if we have one
		while (this->getq (mb) != -1) {
			GadgetContainerMessage<GadgetMessageIdentifier>* mid =
					AsContainerMessage<GadgetMessageIdentifier>(mb);


			if (!mid) {
				ACE_DEBUG ((LM_ERROR, ACE_TEXT ("Invalid message on output queue\n")));
				mb->release();
				return -1;
			}

			//Is this a shutdown message?
			if (mid->getObjectPtr()->id == GADGET_MESSAGE_CLOSE) {
				socket_->send_n(mid->getObjectPtr(),sizeof(GadgetMessageIdentifier));
				return 0;
			}

			GadgetMessageWriter* w = writers_.find(mid->getObjectPtr()->id);

			if (!w) {
				ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unrecognized Message ID received: %d\n"),mid->getObjectPtr()->id));
				mb->release();
				return -1;
			}

			if (w->write(socket_,mb->cont()) < 0) {
				ACE_DEBUG ( (LM_DEBUG, ACE_TEXT ("(%P|%t) Failed to write message to Gadgetron\n")) );
				mb->release ();
				return -1;
			}

			mb->release();
		}

		return 0;

	}


protected:
	ACE_SOCK_Stream* socket_;
	GadgetronSlotContainer<GadgetMessageWriter> writers_;
};

class EXPORTGADGETTOOLS GadgetronConnector: public ACE_Svc_Handler<ACE_SOCK_STREAM, ACE_MT_SYNCH> {

public:
	GadgetronConnector();
	virtual ~GadgetronConnector();

	int open (std::string hostname, std::string port);
	virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE);
	//virtual int handle_output (ACE_HANDLE fd = ACE_INVALID_HANDLE);
	virtual int handle_close (ACE_HANDLE handle, ACE_Reactor_Mask close_mask);
	virtual int svc(void);

	int putq  (  ACE_Message_Block * mb ,  ACE_Time_Value *  timeout = 0) {
		return writer_task_.putq(mb,timeout);
	}

	virtual int process(unsigned int messageid, ACE_Message_Block* mb) {
		mb->release();
		return 0;
	}

	int register_reader(unsigned int slot, GadgetMessageReader* reader);
	int register_writer(unsigned int slot, GadgetMessageWriter* writer) {
		return writer_task_.register_writer(slot,writer);
	}

	int send_gadgetron_configuration_file(std::string config_xml_name);
	int send_gadgetron_configuration_script(std::string config_xml_name);
	int send_gadgetron_parameters(std::string xml_string);

protected:
	//ACE_Reactor_Notification_Strategy notifier_;
	std::string hostname_;
	std::string port_;

	GadgetronSlotContainer<GadgetMessageReader> readers_;
	WriterTask writer_task_;
	//GadgetronSlotContainer<GadgetMessageWriter> writers_;

};

#endif /* GADGETRONCONNECTOR_H_ */
