/** \file gpuoperators_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef GPUOPERATORS_EXPORT_H_
#define GPUOPERATORS_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPU__)
#define EXPORTGPUOPERATORS __declspec(dllexport)
#else
#define EXPORTGPUOPERATORS __declspec(dllimport)
#endif
#else
#define EXPORTGPUOPERATORS
#endif

#endif /* GPUOPERATORS_EXPORT_H_ */
