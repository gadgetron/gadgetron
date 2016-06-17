/** \file   cpusolver_export.h
    \brief  Required definitions for Windows, importing/exporting dll symbols 
    \author Hui Xue
*/

#ifndef CPUSOLVER_EXPORT_H_
#define CPUSOLVER_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CPUSOLVERS__) || defined (cpusolver_EXPORTS)
        #define EXPORTCPUSOLVER __declspec(dllexport)
    #else
        #define EXPORTCPUSOLVER __declspec(dllimport)
    #endif
#else
#define EXPORTCPUSOLVER
#endif

#endif // CPUSOLVER_EXPORT_H_
