/** \file cpufft_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef CPUFFT_EXPORT_H_
#define CPUFFT_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CPUFFT__) || defined (cpufft_EXPORTS)
        #define EXPORTCPUFFT __declspec(dllexport)
    #else
        #define EXPORTCPUFFT __declspec(dllimport)
    #endif
#else
#define EXPORTCPUFFT
#endif

#endif /* CPUCORE_EXPORT_H_ */
