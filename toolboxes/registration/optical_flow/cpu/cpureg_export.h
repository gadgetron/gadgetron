#ifndef _CPUREG_EXPORT_H_
#define _CPUREG_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CPUREG__) || defined (cpureg_EXPORTS)
        #define EXPORTCPUREG __declspec(dllexport)
    #else
        #define EXPORTCPUREG __declspec(dllimport)
    #endif
#else
#define EXPORTCPUREG
#endif

#endif /* _CPUREG_EXPORT_H_ */
