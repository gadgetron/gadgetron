#ifndef MatlabCommunicator_H
#define MatlabCommunicator_H

#include <ace/Synch.h>
#include <ace/Mutex.h>
#include "ace/Event_Handler.h"
//#include "ace/FIFO_Recv_Msg.h"
//#include "ace/SPIPE_Addr.h"
#include "ace/SOCK_Acceptor.h"
#include "ace/Log_Msg.h"
#include "ace/OS_NS_stdio.h"
#include "ace/OS_NS_unistd.h"

#include "Gadget.h"

#include <map>
#include <string>

#include "gadgetronmatlab_export.h"



class EXPORTGADGETSMATLAB MatlabCommunicator : public ACE_Event_Handler
{

public:
	static MatlabCommunicator* instance();

	virtual ACE_HANDLE get_handle (void) const;
	virtual int handle_input (ACE_HANDLE fd);

	//void register_gadget(Gadget* g);
	bool message_gadget(std::string g, ACE_Message_Block* m);

private:
	MatlabCommunicator();
	~MatlabCommunicator();

	static MatlabCommunicator* instance_;
	ACE_Thread_Mutex mutex_;
	ACE_SOCK_Acceptor acceptor_;
};


#endif
