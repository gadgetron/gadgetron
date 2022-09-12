/** \file   pr_export.h
    \brief  Required definitions for Windows, importing/exporting dll symbols 
    \author Hui Xue
*/

#ifndef PR_EXPORT_H_
#define PR_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_PR__) || defined (gadgetron_toolbox_pr_EXPORTS)
        #define EXPORTPR __declspec(dllexport)
    #else
        #define EXPORTPR __declspec(dllimport)
    #endif
#else
    #define EXPORTPR
#endif

#endif // PR_EXPORT_H_
