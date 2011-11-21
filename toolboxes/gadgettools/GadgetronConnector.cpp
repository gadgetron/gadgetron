/*
 * GadgetronConnector.cpp
 *
 *  Created on: Nov 1, 2011
 *      Author: hansenms
 */

#include <ace/SOCK_Connector.h>

#include "GadgetronConnector.h"

#include <iostream>

#define MAXHOSTNAMELENGTH 1024

GadgetronConnector::GadgetronConnector()
	: notifier_ (0, this, ACE_Event_Handler::WRITE_MASK)
{

}

GadgetronConnector::~GadgetronConnector() {
	readers_.clear();
	writers_.clear();
}

int GadgetronConnector::open(std::string hostname, std::string port)
{
	hostname_= hostname;
	port_ = port;

	//Make sure we have a reactor, otherwise assign one from the singleton instance
	if (!this->reactor()) {
		ACE_DEBUG((LM_INFO, ACE_TEXT("Setting reactor")));
		this->reactor(ACE_Reactor::instance());
	}

	//We will add a notification strategy to the message queue to make sure than handle_output gets triggered when packages are on the queue
	this->notifier_.reactor (this->reactor ());
	this->msg_queue ()->notification_strategy (&this->notifier_);

	ACE_INET_Addr server(port_.c_str(),hostname_.c_str());
	ACE_SOCK_Connector connector;

	if (connector.connect(this->peer(),server) == -1) {
		ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"), ACE_TEXT("connect")), -1);
	}

	ACE_TCHAR peer_name[MAXHOSTNAMELENGTH];
	ACE_INET_Addr peer_addr;
	if (peer().get_remote_addr (peer_addr) == 0 && peer_addr.addr_to_string (peer_name, MAXHOSTNAMELENGTH) == 0) {
		ACE_DEBUG ((LM_DEBUG, ACE_TEXT ("(%P|%t) Connection from %s\n"), peer_name));
	}


	if (this->reactor ()->register_handler(this, ACE_Event_Handler::READ_MASK | ACE_Event_Handler::WRITE_MASK) != 0) {
		ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"), ACE_TEXT("Registering read handler")), -2);;
	}

	this->msg_queue ()->notification_strategy (0);

	return this->activate( THR_NEW_LWP | THR_JOINABLE, 1); //Run single threaded. TODO: Add multithreaded support
}


int GadgetronConnector::handle_input(ACE_HANDLE fd)
{
	ssize_t recv_count = 0;
	GadgetMessageIdentifier mid;

	if ((recv_count = peer().recv_n(&mid, sizeof(GadgetMessageIdentifier))) <= 0) {
		ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetronConnector, failed to read message identifier\n")) );
		return -1;
	}

	//Is this a shutdown message?
	if (mid.id == GADGET_MESSAGE_CLOSE) {
		return close();
	}

	GadgetMessageReader* r = readers_.find(mid.id);
	if (r == 0) {
		ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetronConnector, Unknown message id %d received\n"), mid.id) );
		return -1;
	}

	ACE_Message_Block* mb = r->read(&peer());

	if (!mb) {
		ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetronConnector, Failed to read message\n")) );
		return -1;
	}	else {
		if (process(mid.id, mb) < 0) {
			ACE_DEBUG( (LM_ERROR, ACE_TEXT("%P, %l, GadgetronConnector, Failed to process message\n")) );
			return -1;
		}
	}

	return 0;
}

int GadgetronConnector::handle_output(ACE_HANDLE fd)
{
	ACE_Message_Block *mb = 0;
	ACE_Time_Value nowait (ACE_OS::gettimeofday ());

	//Send a package if we have one
	if (-1 != this->getq (mb, &nowait)) {
		GadgetContainerMessage<GadgetMessageIdentifier>* mid =
				AsContainerMessage<GadgetMessageIdentifier>(mb);


		if (!mid) {
			ACE_DEBUG ((LM_ERROR, ACE_TEXT ("Invalid message on output queue\n")));
			mb->release();
			return -1;
		}

		//Is this a shutdown message?
		if (mid->getObjectPtr()->id == GADGET_MESSAGE_CLOSE) {
			peer().send_n(mid->getObjectPtr(),sizeof(GadgetMessageIdentifier));
			return 0;
		}


		GadgetMessageWriter* w = writers_.find(mid->getObjectPtr()->id);

		if (!w) {
			ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unrecognized Message ID received: %d\n"),mid->getObjectPtr()->id));
			return -1;
		}

		if (w->write(&peer(),mb->cont()) < 0) {
			ACE_DEBUG ( (LM_DEBUG, ACE_TEXT ("(%P|%t) Failed to write message to Gadgetron\n")) );
			mb->release ();
			return -1;
		}

		mb->release();
	}

	if (this->msg_queue ()->is_empty ()) {
		//No point in coming back to handle_ouput until something is put on the queue,
		//in which case, the msg queue's notification strategy will tell us

		//Stop the WRITE trigger from the socket
		this->reactor ()->cancel_wakeup(this, ACE_Event_Handler::WRITE_MASK);

		//Get a trigger when stuff is on the queue instead
		this->msg_queue ()->notification_strategy (&this->notifier_);
	} else {
		//There is still more on the queue, let's come back when idle

		//Make sure that we get a wake up when it is possible to write
		this->reactor ()->schedule_wakeup(this, ACE_Event_Handler::WRITE_MASK);

		//Don't wake up from the queue, it may not be possible to write.
		this->msg_queue ()->notification_strategy (0);

		//this->reactor ()->cancel_wakeup(this->notifier_.event_handler(), ACE_Event_Handler::WRITE_MASK);
	}

	return 0;
}

int GadgetronConnector::handle_close(ACE_HANDLE handle, ACE_Reactor_Mask close_mask)
{
	ACE_DEBUG ((LM_INFO, ACE_TEXT ("(%P|%t) Handling close...\n")));
	this->reactor()->end_reactor_event_loop();
	return 0;//this->wait();
}

int GadgetronConnector::svc(void)
{
	//ACE_thread_t old_owner;

	//Take ownership of Reactor
	this->reactor()->owner(ACE_Thread::self ());//, &old_owner);

	this->reactor()->reset_event_loop();

	//Handle the events
	this->reactor()->run_reactor_event_loop();

	//this->reactor()->owner(&old_owner);

	ACE_DEBUG ((LM_DEBUG, ACE_TEXT ("(%P|%t) svc done...\n")));

	return 0;
}

int GadgetronConnector::register_reader(unsigned int slot, GadgetMessageReader *reader)
{
	return readers_.insert(slot,reader);
}


int GadgetronConnector::register_writer(unsigned int slot, GadgetMessageWriter *writer)
{
	return writers_.insert(slot,writer);
}

int GadgetronConnector::send_gadgetron_configuration(std::string config_xml_name)
{
	GadgetMessageIdentifier id;
    id.id = GADGET_MESSAGE_CONFIG_FILE;

    GadgetMessageConfigurationFile ini;
    ACE_OS_String::strncpy(ini.configuration_file, config_xml_name.c_str(),1024);


    if (this->peer().send_n(&id, sizeof(GadgetMessageIdentifier)) != sizeof(GadgetMessageIdentifier)) {
        ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send GadgetMessageIdentifier\n")));
        return -1;
    }

    if (this->peer().send_n(&ini, sizeof(GadgetMessageConfigurationFile)) != sizeof(GadgetMessageConfigurationFile)) {
        ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send GadgetMessageIdentifier\n")));
        return -1;
    }

    return 0;
}



int GadgetronConnector::send_gadgetron_parameters(std::string xml_string)
{
	GadgetMessageIdentifier id;
	id.id = GADGET_MESSAGE_PARAMETER_SCRIPT;

	GadgetMessageScript conf;
	conf.script_length = xml_string.size()+1;
	if (this->peer().send_n(&id, sizeof(GadgetMessageIdentifier)) != sizeof(GadgetMessageIdentifier)) {
		ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send GadgetMessageIdentifier\n")));
		return -1;
	}

	if (this->peer().send_n(&conf, sizeof(GadgetMessageScript)) != sizeof(GadgetMessageScript)) {
		ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send GadgetMessageScript\n")));
		return -1;
	}

	if (this->peer().send_n(xml_string.c_str(), conf.script_length) != conf.script_length) {
		ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send parameter xml\n")));
		return -1;
	}

	return 0;
}









