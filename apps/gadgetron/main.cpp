#include <ace/Log_Msg.h>
#include <ace/Service_Config.h>
#include <ace/Reactor.h>
#include <ace/Get_Opt.h>
#include <ace/OS_NS_string.h>

#include "GadgetServerAcceptor.h"
#include "FileInfo.h"

#include "gadgetron.hxx" //Generated header file for XML configuration
#include <iostream>

#include "url_encode.h"

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

	char * gadgetron_home = ACE_OS::getenv("GADGETRON_HOME");

	if (std::string(gadgetron_home).size() == 0) {
		ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("GADGETRON_HOME variable not set.\n")),-1);
	}

	std::string gcfg = std::string(gadgetron_home) + std::string("/config/gadgetron.xml");

	if (!FileInfo(gcfg).exists()) {
		ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("Gadgetron configuration file %s not found.\n"), gcfg.c_str()),-1);
	}

	ACE_TCHAR schema_file_name[4096];
	ACE_OS::sprintf(schema_file_name, "%s/schema/gadgetron.xsd", gadgetron_home);

	std::string tmp(schema_file_name);
	tmp = url_encode(tmp);
	ACE_OS_String::strncpy(schema_file_name,tmp.c_str(), 4096);

	xml_schema::properties props;
	props.schema_location (
	  "http://gadgetron.sf.net/gadgetron",
	  std::string (schema_file_name));

	ACE_TCHAR port_no[1024];
	try {
		std::auto_ptr<gadgetron::gadgetronConfiguration> cfg(gadgetron::gadgetronConfiguration_(gcfg,0,props));
		ACE_OS_String::strncpy(port_no, cfg->port().c_str(), 1024);
	}  catch (const xml_schema::exception& e) {
		std::cerr << e << std::endl;
		ACE_DEBUG(( LM_DEBUG, ACE_TEXT("XML Parse Error: %s\n"), e.what() ));
		ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("Error parsing configuration file %s.\n"), gcfg.c_str()),-1);
	}

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
