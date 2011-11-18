/*
 * gpucore_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GPUCORE_EXPORT_H_
#define GPUCORE_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPUCORE__) || defined (gpucore_EXPORTS)
#define EXPORTGPUCORE __declspec(dllexport)
#else
#define EXPORTGPUCORE __declspec(dllimport)
#endif
#else
#define EXPORTGPUCORE
#endif

#endif /* GPUCORE_EXPORT_H_ */
