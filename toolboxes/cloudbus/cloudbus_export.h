#ifndef CLOUDBUS_EXPORT_H_
#define CLOUDBUS_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CLOUDBUS__) || defined (gadgetron_toolbox_cloudbus_EXPORTS)
        #define EXPORTCLOUDBUS __declspec(dllexport)
    #else
        #define EXPORTCLOUDBUS __declspec(dllimport)
    #endif
#else
    #define EXPORTCLOUDBUS
#endif

#endif
