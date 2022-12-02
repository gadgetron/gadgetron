/** \file   cpuwavelet_export.h
    \brief  Required definitions for Windows, importing/exporting dll symbols 
    \author Hui Xue
*/

#ifndef CPUWAVELET_EXPORT_H_
#define CPUWAVELET_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CPUDWT__) || defined (cpudwt_EXPORTS)
        #define EXPORTCPUDWT __declspec(dllexport)
    #else
        #define EXPORTCPUDWT __declspec(dllimport)
    #endif
#else
#define EXPORTCPUDWT
#endif

#endif // CPUWAVELET_EXPORT_H_
