#include "GadgetServerAcceptor.h"
#include "FileInfo.h"
#include "url_encode.h"
#include "gadgetron_xml.h"
#include "gadgetron_config.h"
#include "gadgetron_paths.h"
#include "CloudBus.h"
#include "gadgetron_rest.h"
#include "gadgetron_system_info.h"

#include <ace/Log_Msg.h>
#include <ace/Service_Config.h>
#include <ace/Reactor.h>
#include <ace/Get_Opt.h>
#include <ace/OS_NS_string.h>
#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>


#ifdef _WIN32
#include <windows.h>
#include <Shlwapi.h>
#pragma comment(lib, "shlwapi.lib")
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif // _WIN32

#include <boost/filesystem.hpp>
using namespace boost::filesystem;

using namespace Gadgetron;

#define GT_WORKING_DIRECTORY "workingDirectory"

namespace Gadgetron {


  bool create_folder_with_all_permissions(const std::string& workingdirectory)
  {
    if ( !boost::filesystem::exists(workingdirectory) )
      {
        boost::filesystem::path workingPath(workingdirectory);
        try
          {
            boost::filesystem::create_directories(workingPath);
          }
        catch (...)
	  {
	    GERROR("Error creating the working directory.\n");
	    return false;
	  }

        // set the permission for the folder
#ifdef _WIN32
	try
	  {
	    boost::filesystem::permissions(workingPath, all_all);
	  }
	catch(...)
	  {
	    GERROR("Error changing the permission of the working directory.\n");
	    return false;
	  }
#else
	// in case an older version of boost is used in non-win system
	// the system call is used
	int res = chmod(workingPath.string().c_str(), S_IRUSR|S_IWUSR|S_IXUSR|S_IRGRP|S_IWGRP|S_IXGRP|S_IROTH|S_IWOTH|S_IXOTH);
	if ( res != 0 )
	  {
	    GERROR("Error changing the permission of the working directory.\n");
	    return false;
	  }
#endif // _WIN32
      }

    return true;
  }

}

void print_usage()
{
  GINFO("Usage: \n");
  GINFO("gadgetron   -p <PORT>                      (default 9002)       \n");
  GINFO("            -r <RELAY HOST>                (default localhost)  \n");
  GINFO("            -l <RELAY PORT>                (default 0, disabled)\n");
}

int ACE_TMAIN(int argc, ACE_TCHAR *argv[])
{
  std::string  gadgetron_home = get_gadgetron_home();

  if (gadgetron_home.size() == 0) {
    GERROR("GADGETRON_HOME variable not set.\n");
    return -1;
  }

  std::string gcfg = gadgetron_home + std::string("/") + std::string(GADGETRON_CONFIG_PATH) + std::string("/gadgetron.xml");
  if (!FileInfo(gcfg).exists()) {
    GERROR("Gadgetron configuration file %s not found.\n", gcfg.c_str());
    return -1;
  }


  ACE_TCHAR port_no[1024];
  ACE_TCHAR relay_host[1024];
  uint16_t  relay_port = 0;

  ACE_OS_String::strncpy(relay_host, "localhost", 1024);

  std::map<std::string, std::string> gadget_parameters;

  // the working directory of gadgetron should always be set
  bool workingDirectorySet = false;

  GadgetronXML::GadgetronConfiguration c;
  try
    {
      std::ifstream t(gcfg.c_str());
      std::string gcfg_text((std::istreambuf_iterator<char>(t)),
			    std::istreambuf_iterator<char>());
      
      GadgetronXML::deserialize(gcfg_text.c_str(), c);
      ACE_OS_String::strncpy(port_no, c.port.c_str(), 1024);

      if (c.cloudBus) {
	ACE_OS_String::strncpy(relay_host, c.cloudBus->relayAddress.c_str(), 1024);
	relay_port = c.cloudBus->port;
      }

      if (c.rest) {
	Gadgetron::ReST::GadgetronReST::port_ = c.rest->port;
	Gadgetron::ReST::GadgetronReST::instance(); //Fire up the server
	
	Gadgetron::ReST::GadgetronReST::instance()->resource["^/info$"]["GET"]=[](Gadgetron::ReST::GadgetronReST::Response& response,
										  std::shared_ptr<Gadgetron::ReST::GadgetronReST::Request> request)
	  {
	    std::stringstream ss;
	    print_system_information(ss);
	    std::string content = ss.str();
	    response << "HTTP/1.1 200 OK\r\nContent-Length: " << content.length() << "\r\n\r\n" << content;
	  };
	Gadgetron::ReST::GadgetronReST::instance()->update_resources();
      }
      
      for (std::vector<GadgetronXML::GadgetronParameter>::iterator it = c.globalGadgetParameter.begin();
	   it != c.globalGadgetParameter.end();
	   ++it)
	{
	  std::string key = it->name;
	  std::string value = it->value;
      
	  gadget_parameters[key] = value;
	  
	  if ( key == std::string(GT_WORKING_DIRECTORY) ) workingDirectorySet = true;
        }
    }  catch (std::runtime_error& e) {
    GERROR("XML Parse Error: %s\n", e.what());
    GERROR("Error parsing configuration file %s.\n", gcfg.c_str());
    return -1;
  }

  static const ACE_TCHAR options[] = ACE_TEXT(":p:r:l:");
  ACE_Get_Opt cmd_opts(argc, argv, options);

  int option;
  while ((option = cmd_opts()) != EOF) {
    switch (option) {
    case 'p':
      ACE_OS_String::strncpy(port_no, cmd_opts.opt_arg(), 1024);
      break;
    case 'r':
      ACE_OS_String::strncpy(relay_host, cmd_opts.opt_arg(), 1024);
      break;
    case 'l':
      relay_port = std::atoi(cmd_opts.opt_arg());
      break;
    case ':':
      print_usage();
      GERROR("-%c requires an argument.\n", cmd_opts.opt_opt());
      return -1;
      break;
    default:
      print_usage();
      GERROR("Command line parse error\n");
      return -1;
      break;
    }
  }

  if (relay_port > 0) {
    GINFO("Starting cloudBus: %s:%d\n", relay_host, relay_port);
    Gadgetron::CloudBus::set_relay_address(relay_host);
    Gadgetron::CloudBus::set_relay_port(relay_port);
    Gadgetron::CloudBus::set_gadgetron_port(std::atoi(port_no));
    Gadgetron::CloudBus* cb = Gadgetron::CloudBus::instance();//This actually starts the bus.
    gadget_parameters["using_cloudbus"] = std::string("true"); //This is our message to the Gadgets that we have activated the bus
  }


  // if the working directory is not set, use the default path
  if ( !workingDirectorySet )
    {
#ifdef _WIN32
      gadget_parameters[std::string(GT_WORKING_DIRECTORY)] = std::string("c:\\temp\\gadgetron\\");
#else
      gadget_parameters[std::string(GT_WORKING_DIRECTORY)] = std::string("/tmp/gadgetron/");
#endif // _WIN32
    }

  // check and create workingdirectory
  std::string workingDirectory = gadget_parameters[std::string(GT_WORKING_DIRECTORY)];
  if ( !Gadgetron::create_folder_with_all_permissions(workingDirectory) )
    {
      GERROR("Gadgetron creating working directory %s failed ... \n", workingDirectory.c_str());
      return -1;
    }

  GINFO("Configuring services, Running on port %s\n", port_no);

  ACE_INET_Addr port_to_listen (port_no);
  GadgetServerAcceptor acceptor;
  acceptor.global_gadget_parameters_ = gadget_parameters;
  acceptor.reactor (ACE_Reactor::instance ());
  if (acceptor.open (port_to_listen) == -1)
    return 1;

  ACE_Reactor::instance()->run_reactor_event_loop ();

  return 0;
}
