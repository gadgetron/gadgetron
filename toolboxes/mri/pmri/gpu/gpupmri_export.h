/** \file gpupmri_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef GPUPMRI_EXPORT_H_
#define GPUPMRI_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPU__) || defined (gpuparallelmri_EXPORTS)
#define EXPORTGPUPMRI __declspec(dllexport)
#else
#define EXPORTGPUPMRI __declspec(dllimport)
#endif
#else
#define EXPORTGPUPMRI
#endif


#endif /* GPUPMRI_EXPORT_H_ */
