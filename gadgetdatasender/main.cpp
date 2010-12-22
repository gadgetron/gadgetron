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
  static const ACE_TCHAR options[] = ACE_TEXT(":p:h:c:f:");
  
  ACE_Get_Opt cmd_opts(argc, argv, options);
  
  ACE_TCHAR port_no[1024];
  ACE_OS_String::strncpy(port_no, "9002", 1024);
  
  ACE_TCHAR hostname[1024];
  ACE_OS_String::strncpy(hostname, "localhost", 1024);
  
  ACE_TCHAR filename[4096];
  ACE_OS_String::strncpy(filename, "./data.dat", 4096);

  ACE_TCHAR config_lib[1024];
  ACE_OS_String::strncpy(config_lib, "core", 1024);

  ACE_TCHAR config_name[1024];
  ACE_OS_String::strncpy(config_name, "default", 1024);
 
   int option;
  while ((option = cmd_opts()) != EOF) {
    switch (option) {
    case 'p':
      ACE_OS_String::strncpy(port_no, cmd_opts.opt_arg(), 1024);
      break;
    case 'h':
      ACE_OS_String::strncpy(hostname, cmd_opts.opt_arg(), 1024);
      break;
    case 'f':
      ACE_OS_String::strncpy(filename, cmd_opts.opt_arg(), 4096);
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
  ACE_OS_String::strncpy(conf.configurator_lib,config_lib,1024);
  ACE_OS_String::strncpy(conf.configurator_name,config_name,1024);
  conf.configuration_length = 1024;

  //Let's declare an empty config file for now.
  char configuration_file[1024];
  ACE_OS::memset(configuration_file,0,1024);

  ACE_DEBUG(( LM_INFO, ACE_TEXT("Sending configuration %s@%s \n"), conf.configurator_name, conf.configurator_lib ));
  
  ACE_INET_Addr server(port_no,hostname);
  ACE_SOCK_Connector connector;
  ACE_SOCK_Stream peer;
  
  if (connector.connect(peer,server) == -1) {
    ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"), ACE_TEXT("connect")), 100);
  } else {
    peer.send_n(&id, sizeof(GadgetMessageIdentifier));
    peer.send_n(&conf, sizeof(GadgetMessageConfigurator));
    peer.send_n(configuration_file, conf.configuration_length);
    peer.close();
  }
  
  
  return 0;
}
