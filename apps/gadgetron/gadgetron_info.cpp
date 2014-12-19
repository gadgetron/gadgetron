//#include "ace/OS_NS_stdlib.h"
#include "ace/OS_NS_string.h"
//#include "ace/OS_NS_stdio.h"
#include "ace/DLL.h"
#include "ace/DLL_Manager.h"
//#include "ace/OS_NS_netdb.h"

#include "gadgetron_config.h"
#include "Gadget.h"

#include <iostream>


using namespace Gadgetron;

int main(int argc, char** argv)
{
  std::cout << "Gadgetron Version Info" << std::endl;
  std::cout << "  -- Version  : " << GADGETRON_VERSION_STRING << std::endl;
  std::cout << "  -- Git SHA1 : " << GADGETRON_GIT_SHA1_HASH << std::endl << std::endl;
  
  if (argc == 1) {
    return 0;
  }

  if ((argc == 2) || (argc > 3)) {
    std::cout << "Invalid number of arguments (" << argc -1 << ")." << std::endl;
    std::cout << "Usage (gadget library info):  " << argc << std::endl;
    std::cout << " -- gadgetron_info <SHARED LIB> <GADGET_INFO>" << std::endl;
    return -1; 
  }

  const char* DLL = argv[1];
  const char* component_name = argv[2];

  //We must be investigating a certain gadget
  std::cout << "Examining Gadget (SHARED LIB): " << component_name << " (" << DLL << ")" << std::endl;

  //Attempt to load Gadget
  //ACE_DLL_Manager* dllmgr = ACE_DLL_Manager::instance();
  
  ACE_DLL_Handle dll;// = 0;
  ACE_SHLIB_HANDLE dll_handle = 0;
  
  ACE_TCHAR dllname[1024];
#if defined(WIN32) && defined(_DEBUG)
  ACE_OS::sprintf(dllname, "%s%sd",ACE_DLL_PREFIX, DLL);
#else
  ACE_OS::sprintf(dllname, "%s%s",ACE_DLL_PREFIX, DLL);
#endif

  ACE_TCHAR factoryname[1024];
  ACE_OS::sprintf(factoryname, "make_%s", component_name);
  
  if (dll.open(dllname, ACE_DEFAULT_SHLIB_MODE, dll_handle )) {
    std::cout << "Failed to load DLL (" << DLL << "), Possible reasons:" << std::endl;
    std::cout << "   - Name of DLL is wrong" << std::endl;
    std::cout << "   - Path of DLL is not in your DLL search path (LD_LIBRARY_PATH on Unix)" << std::endl;
    std::cout << "   - Path of other DLLs that this DLL depends on is not in the search path" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Set environment variable ACE_DEBUG=1 to get more information" << std::endl << std::endl; 
    return 0;
  } 

  //Function pointer
  typedef Gadget* (*ComponentCreator) (void);

  void *void_ptr = dll.symbol (factoryname);
  ptrdiff_t tmp = reinterpret_cast<ptrdiff_t> (void_ptr);
  ComponentCreator cc = reinterpret_cast<ComponentCreator> (tmp);
  
  if (cc == 0) {
    std::cout << "Failed to load factory (" << factoryname << ") from DLL (" << dllname << ")" << std::endl;
    return -1;
  }
  
  Gadget* g = cc();
  if (!g) {
    std::cout << "Failed to create component using factory" << std::endl;
    return 0;
  }

  std::cout << "  -- Gadget compiled against Gadgetron version " << g->get_gadgetron_version() << std::endl;

  delete g;

  return 0;
}
