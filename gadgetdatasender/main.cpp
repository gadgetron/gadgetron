#include "ace/INET_Addr.h"
#include "ace/SOCK_Stream.h"
#include "ace/SOCK_Connector.h"
#include "ace/Log_Msg.h"
#include "ace/Get_Opt.h"
#include "ace/OS_String.h"

#include "../gadgetheaders.h"
#include "siemensraw.hpp"

int ACE_TMAIN(int argc, ACE_TCHAR *argv[] )
{
  static const ACE_TCHAR options[] = ACE_TEXT(":p:h:c:");
  
  ACE_Get_Opt cmd_opts(argc, argv, options);
  
  ACE_TCHAR port_no[1024];
  ACE_OS_String::strncpy(port_no, "9002", 1024);
  
  ACE_TCHAR hostname[1024];
  ACE_OS_String::strncpy(hostname, "localhost", 1024);
  
  int configuration = 0;

  int option;
  while ((option = cmd_opts()) != EOF) {
    switch (option) {
    case 'p':
      ACE_OS_String::strncpy(port_no, cmd_opts.opt_arg(), 1024);
      break;
    case 'h':
      ACE_OS_String::strncpy(hostname, cmd_opts.opt_arg(), 1024);
      break;
    case 'c':
      configuration = ACE_OS::atoi(cmd_opts.opt_arg());
      break;
    case ':':
      ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("-%c requires an argument.\n"), cmd_opts.opt_opt()),-1);
      break;
    default:
      ACE_ERROR_RETURN( (LM_ERROR, ACE_TEXT("Command line parse error\n")), -1);
    }    
  }
  
  ACE_DEBUG(( LM_INFO, ACE_TEXT("Gadgetron Data Sender\n") ));
  
  GadgetMessageIdentifier id;
  GadgetMessageConfigurator conf;

  id.id = GADGET_MESSAGE_CONFIGURATION;

  if (configuration > 0) {
    ACE_DEBUG(( LM_INFO, ACE_TEXT("Sending configuration %d \n"), configuration ));
    
    ACE_INET_Addr server(port_no,hostname);
    ACE_SOCK_Connector connector;
    ACE_SOCK_Stream peer;
    
    if (connector.connect(peer,server) == -1) {
      ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"), ACE_TEXT("connect")), 100);
    } else {
      peer.send_n(&id, sizeof(GadgetMessageIdentifier));
      peer.close();
    }
    
  } else {	
    ACE_DEBUG(( LM_INFO, ACE_TEXT("Configuration %d is not valid\n"), configuration ));
    
    
  }
  
  return 0;
}
