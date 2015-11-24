/** \file   cpuwavelet_export.h
    \brief  Required definitions for Windows, importing/exporting dll symbols 
    \author Hui Xue
*/

#ifndef CPUWAVELET_EXPORT_H_
#define CPUWAVELET_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CPUWAVELET__) || defined (cpuwavelet_EXPORTS)
        #define EXPORTCPUWAVELET __declspec(dllexport)
    #else
        #define EXPORTCPUWAVELET __declspec(dllimport)
    #endif
#else
#define EXPORTCPUWAVELET
#endif

#endif // CPUWAVELET_EXPORT_H_
