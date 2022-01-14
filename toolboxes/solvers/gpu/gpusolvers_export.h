/** \file gpusolvers_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef GPUSOLVERS_EXPORT_H_
#define GPUSOLVERS_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPU__) || defined (gpusolvers_EXPORTS)
#define EXPORTGPUSOLVERS __declspec(dllexport)
#else
#define EXPORTGPUSOLVERS __declspec(dllimport)
#endif
#else
#define EXPORTGPUSOLVERS
#endif

#endif /* GPUSOLVERS_EXPORT_H_ */
