#include <ace/Log_Msg.h>
#include <ace/Service_Config.h>
#include <ace/Reactor.h>
#include <ace/Get_Opt.h>
#include <ace/OS_NS_string.h>

#include "GadgetServerAcceptor.h"

void print_usage()
{
	ACE_DEBUG((LM_INFO, ACE_TEXT("Usage: \n") ));
	ACE_DEBUG((LM_INFO, ACE_TEXT("gadgetron   -p <PORT>                      (default 9002)       \n") ));
}

int ACE_TMAIN(int argc, ACE_TCHAR *argv[])
{
	ACE_TRACE(( ACE_TEXT("main") ));

	ACE_LOG_MSG->priority_mask( LM_INFO | LM_NOTICE | LM_ERROR| LM_DEBUG,
			ACE_Log_Msg::PROCESS);

	ACE_TCHAR port_no[1024];
	ACE_OS_String::strncpy(port_no, "9002", 1024);

	static const ACE_TCHAR options[] = ACE_TEXT(":p:");
	ACE_Get_Opt cmd_opts(argc, argv, options);

	int option;
	while ((option = cmd_opts()) != EOF) {
		switch (option) {
		case 'p':
			ACE_OS_String::strncpy(port_no, cmd_opts.opt_arg(), 1024);
			break;
		case ':':
			print_usage();
			ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("-%c requires an argument.\n"), cmd_opts.opt_opt()),-1);
			break;
		default:
			print_usage();
			ACE_ERROR_RETURN( (LM_ERROR, ACE_TEXT("Command line parse error\n")), -1);
			break;
		}
	}


	ACE_DEBUG(( LM_DEBUG, ACE_TEXT("%IConfiguring services, Running on port %s\n"), port_no ));

	ACE_INET_Addr port_to_listen (port_no);
	GadgetServerAcceptor acceptor;
	acceptor.reactor (ACE_Reactor::instance ());
	if (acceptor.open (port_to_listen) == -1)
		return 1;

	ACE_Reactor::instance()->run_reactor_event_loop ();

	return 0;
}
