/*
 * gpupmri_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GPUPMRI_EXPORT_H_
#define GPUPMRI_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPUPMRI__) || defined (gpuparallelmri_EXPORTS)
#define EXPORTGPUPMRI __declspec(dllexport)
#else
#define EXPORTGPUPMRI __declspec(dllimport)
#endif
#else
#define EXPORTGPUPMRI
#endif


#endif /* GPUPMRI_EXPORT_H_ */
