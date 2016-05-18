/** \file   cmr_export.h
    \brief  Required definitions for Windows, importing/exporting dll symbols 
    \author Hui Xue
*/

#ifndef CMR_EXPORT_H_
#define CMR_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CMR__) || defined (gadgetron_toolbox_cmr_EXPORTS)
        #define EXPORTCMR __declspec(dllexport)
    #else
        #define EXPORTCMR __declspec(dllimport)
    #endif
#else
    #define EXPORTCMR
#endif

#endif // CMR_EXPORT_H_
