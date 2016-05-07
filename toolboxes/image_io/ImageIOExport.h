/** \file       ImageIOExport.h
    \brief      Implement export/import for image_io toolbox
    \author     Hui Xue
*/

#pragma once

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_IMAGE_ANALYZE_IO__) || defined (image_analyze_io_EXPORTS)
        #define EXPORTIMAGEIO __declspec(dllexport)
    #else
        #define EXPORTIMAGEIO __declspec(dllimport)
    #endif
#else
    #define EXPORTIMAGEIO
#endif
