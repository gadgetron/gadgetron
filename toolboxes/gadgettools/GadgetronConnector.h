/*
 * GadgetronConnector.h
 *
 *  Created on: Nov 1, 2011
 *      Author: hansenms
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

class EXPORTGADGETTOOLS GadgetronConnector: public ACE_Svc_Handler<ACE_SOCK_STREAM, ACE_MT_SYNCH> {

public:
	GadgetronConnector();
	virtual ~GadgetronConnector();

	int open (std::string hostname, std::string port);
	virtual int handle_input (ACE_HANDLE fd = ACE_INVALID_HANDLE);
	virtual int handle_output (ACE_HANDLE fd = ACE_INVALID_HANDLE);
	virtual int handle_close (ACE_HANDLE handle, ACE_Reactor_Mask close_mask);
	virtual int svc(void);

	virtual int process(unsigned int messageid, ACE_Message_Block* mb) {
			mb->release();
			return 0;
	}

	int register_reader(unsigned int slot, GadgetMessageReader* reader);
	int register_writer(unsigned int slot, GadgetMessageWriter* writer);

	int send_gadgetron_configuration(std::string config_xml_name);
	int send_gadgetron_parameters(std::string xml_string);

protected:
	ACE_Reactor_Notification_Strategy notifier_;
	std::string hostname_;
	std::string port_;

	GadgetronSlotContainer<GadgetMessageReader> readers_;
	GadgetronSlotContainer<GadgetMessageWriter> writers_;

};

#endif /* GADGETRONCONNECTOR_H_ */
