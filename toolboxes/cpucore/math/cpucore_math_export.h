/*
 * ndarray_export.h
 *
 *  Created on: Nov 29, 2011
 *      Author: hansenms
 */

#ifndef CPUCORE_EXPORT_H_
#define CPUCORE_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_CPUCORE__) || defined (cpucore_EXPORTS)
#define EXPORTCPUCORE __declspec(dllexport)
#else
#define EXPORTCPUCORE __declspec(dllimport)
#endif
#else
#define EXPORTCPUCORE
#endif

#endif /* CPUCORE_EXPORT_H_ */
