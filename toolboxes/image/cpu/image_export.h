/** \file   image_export.h
\brief  Required definitions for Windows, importing/exporting dll symbols
\author Hui Xue
*/

#ifndef  IMAGE_EXPORT_H_
#define IMAGE_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_IMAGE__) || defined (image_EXPORTS)
        #define EXPORTIMAGE __declspec(dllexport)
    #else
        #define EXPORTIMAGE __declspec(dllimport)
    #endif
#else
    #define EXPORTIMAGE
#endif

#endif // IMAGE_EXPORT_H_
