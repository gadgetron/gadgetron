#include "GadgetronConnector.h"
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"

#include <ace/SOCK_Acceptor.h>
#include <ace/Addr.h>
#include <ace/INET_Addr.h>
#include <ace/Log_Msg.h>
#include <ace/Get_Opt.h>
#include <ace/OS_NS_string.h>

#include <fstream>
#include <time.h>
#include <iomanip>
#include <iostream>

using namespace Gadgetron;

int ACE_TMAIN(int argc, ACE_TCHAR *argv[] )
{
	GadgetronConnector con;

	std::string host("localhost");
	std::string port("9002");

	ACE_TCHAR hostname[1024];
	//We will do a little trick to figure out what the hostname would be according to ACE
	ACE_SOCK_Acceptor listener (ACE_Addr::sap_any);
	ACE_INET_Addr addr;
	listener.get_local_addr (addr);
	ACE_OS_String::strncpy(hostname, addr.get_host_name(), 1024);

	host = std::string(hostname);
        
	if (argc > 1) {
		host = std::string(argv[1]);
	}

	if (argc > 2) {
		port = std::string(argv[2]);
	}

	if (con.open(host,port) != 0) {
	  GERROR("Unable to connect to the Gadgetron host\n");
	  return -1;
	}

	//Tell Gadgetron which XML configuration to run.
	if (con.send_gadgetron_configuration_file(std::string("isalive.xml")) != 0) {
	  GERROR("Unable to send XML configuration to the Gadgetron host\n");
	  return -1;
	}


	GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
			new GadgetContainerMessage<GadgetMessageIdentifier>();

	m1->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;

	if (con.putq(m1) == -1) {
	  GERROR("Unable to put CLOSE package on queue\n");
	  return -1;
	}

	con.wait();

	return 0;
}
