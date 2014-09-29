#include "GadgetServerAcceptor.h"
#include "FileInfo.h"
#include "url_encode.h"
#include "gadgetron_xml.h"
#include "gadgetron_config.h"
#include "gadgetron_paths.h"

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
        if ( !boost::filesystem::create_directory(workingPath) )
        {
            ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("Error creating the working directory.\n")), false);
        }

        // set the permission for the folder
        #ifdef _WIN32
            try
            {
                boost::filesystem::permissions(workingPath, all_all);
            }
            catch(...)
            {
                ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("Error changing the permission of the working directory.\n")), false);
            }
        #else
            // in case an older version of boost is used in non-win system
            // the system call is used
            int res = chmod(workingPath.string().c_str(), S_IRUSR|S_IWUSR|S_IXUSR|S_IRGRP|S_IWGRP|S_IXGRP|S_IROTH|S_IWOTH|S_IXOTH);
            if ( res != 0 )
            {
                ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("Error changing the permission of the working directory.\n")), false);
            }
        #endif // _WIN32
    }

    return true;
}

}

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

    std::string  gadgetron_home = get_gadgetron_home();

    if (gadgetron_home.size() == 0) {
        ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("GADGETRON_HOME variable not set.\n")),-1);
    }

    std::string gcfg = gadgetron_home + std::string("/") + std::string(GADGETRON_CONFIG_PATH) + std::string("/gadgetron.xml");
    if (!FileInfo(gcfg).exists()) {
        ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("Gadgetron configuration file %s not found.\n"), gcfg.c_str()),-1);
    }


    ACE_TCHAR port_no[1024];
    std::map<std::string, std::string> gadget_parameters;

    // the working directory of gadgetron should always be set
    bool workingDirectorySet = false;

    try
    {
      std::ifstream t(gcfg.c_str());
      std::string gcfg_text((std::istreambuf_iterator<char>(t)),
			    std::istreambuf_iterator<char>());
      
      GadgetronXML::GadgetronConfiguration c;
      GadgetronXML::deserialize(gcfg_text.c_str(), c);
      ACE_OS_String::strncpy(port_no, c.port.c_str(), 1024);
 
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
        ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("Gadgetron creating working directory %s failed ... \n"), workingDirectory.c_str()),-1);
    }

    ACE_DEBUG(( LM_DEBUG, ACE_TEXT("%IConfiguring services, Running on port %s\n"), port_no ));

    ACE_INET_Addr port_to_listen (port_no);
    GadgetServerAcceptor acceptor;
    acceptor.global_gadget_parameters_ = gadget_parameters;
    acceptor.reactor (ACE_Reactor::instance ());
    if (acceptor.open (port_to_listen) == -1)
        return 1;

    ACE_Reactor::instance()->run_reactor_event_loop ();

    return 0;
}
