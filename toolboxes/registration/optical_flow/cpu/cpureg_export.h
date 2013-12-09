#ifndef _CPUREG_EXPORT_H_
#define _CPUREG_EXPORT_H_

#if defined (WIN32)
    #ifdef BUILD_TOOLBOX_STATIC
        #define EXPORTCPUREG
    #else
        #if defined (__BUILD_GADGETRON_CPUREG__) || defined (cpureg_EXPORTS)
            #define EXPORTCPUREG __declspec(dllexport)
        #else
            #define EXPORTCPUREG __declspec(dllimport)
        #endif
    #endif
#else
#define EXPORTCPUREG
#endif

#endif /* _CPUREG_EXPORT_H_ */
