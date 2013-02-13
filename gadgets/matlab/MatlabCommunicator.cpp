#include "MatlabCommunicator.h"

#include <iostream>

MatlabCommunicator* MatlabCommunicator::instance()
{
	if (!instance_) instance_ = new MatlabCommunicator();
	return instance_;
}

MatlabCommunicator::MatlabCommunicator()
: mutex_("MatlabCommunicatorMutex")
{

	std::cout << "Starting up Matlab Communicator" << std::endl;
	const char* gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");
	std::string path_name = std::string(gadgetron_home) + std::string("/matlab");


	const ACE_INET_Addr listen_addr("9003");

	std::cout << "About to open acceptor:" << std::endl;
	if (this->acceptor_.open (listen_addr) == -1)
		ACE_ERROR ((LM_ERROR, ACE_TEXT ("%p\n"), ACE_TEXT ("open")));

	// Register with the Reactor.
	if (ACE_Reactor::instance ()->register_handler
			(this, ACE_Event_Handler::ACCEPT_MASK) == -1)
		ACE_ERROR ((LM_ERROR, ACE_TEXT ("%p\n"), ACE_TEXT ("register_handler")));

}

MatlabCommunicator::~MatlabCommunicator()
{
	this->acceptor_.close ();
}

ACE_HANDLE MatlabCommunicator::get_handle (void) const
{
	//return this->fifo_reader_.get_handle ();
	return this->acceptor_.get_handle ();
}


int MatlabCommunicator::handle_input (ACE_HANDLE)
{

	std::cout << "BLAH BLAAHHH" << std::endl;
	char buf[BUFSIZ];

	ACE_DEBUG ((LM_DEBUG, ACE_TEXT ("handle_input\n")));

	ACE_Str_Buf msg (buf, 0, sizeof buf);

	ACE_SOCK_Stream new_stream;
	if (acceptor_.accept (new_stream, 0) == -1)
		ACE_ERROR_RETURN ((LM_ERROR,
				ACE_TEXT ("%p\n"),
				ACE_TEXT ("accept")),
				1);

	ACE_DEBUG ((LM_DEBUG,
			ACE_TEXT ("Accepted connection\n")));

	int  n;
	while ((n = new_stream.recv (buf, sizeof buf)) > 0)
	{
		ACE_OS::fprintf (stderr,
				"%s\n",
				buf);
		ACE_OS::write (ACE_STDOUT,
				buf,
				n);
	}

	if (n == -1)
	{
		ACE_DEBUG ((LM_DEBUG,
				ACE_TEXT ("End of connection. Closing handle\n")));
		new_stream.close ();
	}

	return 0;
}

bool MatlabCommunicator::message_gadget(std::string gadget, ACE_Message_Block* m)
{
	return false;
}

MatlabCommunicator* MatlabCommunicator::instance_ = NULL;
