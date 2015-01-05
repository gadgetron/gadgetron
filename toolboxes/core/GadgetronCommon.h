#ifndef GADGETRONCOMMON_H
#define GADGETRONCOMMON_H
//#include "log.h"
#include <iostream>

#ifndef _WIN32

#define GCC_VERSION (__GNUC__ * 10000           \
                     + __GNUC_MINOR__ * 1000    \
                     + __GNUC_PATCHLEVEL__)

#if GCC_VERSION < 42000
#pragma message ("GCC version is older than 4.2.0")
#define GCC_OLD_FLAG 1
#endif

#else

    // disable warning 4251, needs to have dll-interface to be used by clients
    #pragma warning( disable : 4251 )

    // warning C4344: behavior change: use of explicit template arguments
    #pragma warning( disable : 4344 )

    // The POSIX name for this item is deprecated. Instead, use the ISO C++ conformant name
    #pragma warning( disable : 4996 )

#endif // _WIN32

#endif  //GADGETRONCOMMON_H
