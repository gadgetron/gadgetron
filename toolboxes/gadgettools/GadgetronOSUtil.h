#ifndef GADGETRONOSUTIL_H_
#define GADGETRONOSUTIL_H_

#include "gadgettools_export.h"
#include "GadgetronCommon.h"
#include <string>
#include <iostream>

#ifdef _WIN32
    #include <windows.h>
#else
    
#endif // _WIN32

namespace Gadgetron{

    EXPORTGADGETTOOLS bool create_folder_with_all_permissions(const std::string& workingdirectory);

}

#endif /* GADGETRONOSUTIL_H_ */
