
#include "GadgetronOSUtil.h"

#ifdef _WIN32
    #include <windows.h>
#else
    #include <sys/types.h>
    #include <sys/stat.h>
#endif // _WIN32

#include <boost/filesystem.hpp>
using namespace boost::filesystem;

namespace Gadgetron{

    bool create_folder_with_all_permissions(const std::string& workingdirectory)
    {
        if ( !boost::filesystem::exists(workingdirectory) )
        {
            boost::filesystem::path workingPath(workingdirectory);
            if ( !boost::filesystem::create_directory(workingPath) )
            {
                GERROR_STREAM("Error creating the working directory " << workingdirectory);
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
                    GERROR_STREAM("Error changing the permission of the working directory " << workingdirectory);
                }
            #else
                // in case an older version of boost is used in non-win system
                // the system call is used
                int res = chmod(workingPath.c_str(), S_IRUSR|S_IWUSR|S_IXUSR|S_IRGRP|S_IWGRP|S_IXGRP|S_IROTH|S_IWOTH|S_IXOTH);
                if ( res != 0 )
                {
                    GERROR_STREAM("Error changing the permission of the working directory " << workingdirectory);
                }
            #endif // _WIN32
        }

        return true;
    }
}
