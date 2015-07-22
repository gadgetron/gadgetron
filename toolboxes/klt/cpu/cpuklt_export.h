/** \file   cpuklt_export.h
    \brief  Required definitions for Windows, importing/exporting dll symbols 
    \author Hui Xue
*/

#ifndef CPUKLT_EXPORT_H_
#define CPUKLT_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CPUKLT__) || defined (cpuklt_EXPORTS)
        #define EXPORTCPUKLT __declspec(dllexport)
    #else
        #define EXPORTCPUKLT __declspec(dllimport)
    #endif
#else
    #define EXPORTCPUKLT
#endif

#endif // CPUKLT_EXPORT_H_
