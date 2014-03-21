/** \file gpufft_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef GPUFFT_EXPORT_H_
#define GPUFFT_EXPORT_H_

#if defined (WIN32)
    #ifdef BUILD_TOOLBOX_STATIC
        #define EXPORTGPUFFT
    #else
        #if defined (__BUILD_GADGETRON_GPUFFT__) || defined (gpufft_EXPORTS)
            #define EXPORTGPUFFT __declspec(dllexport)
        #else
            #define EXPORTGPUFFT __declspec(dllimport)
        #endif
    #endif
#else
#define EXPORTGPUFFT
#endif

#endif
