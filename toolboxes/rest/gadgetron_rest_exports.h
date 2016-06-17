#ifndef REST_EXPORT_H_
#define REST_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_REST__) || defined (gadgetron_toolbox_rest_EXPORTS)
        #define EXPORTREST __declspec(dllexport)
    #else
        #define EXPORTREST __declspec(dllimport)
    #endif
#else
    #define EXPORTREST
#endif

#endif
