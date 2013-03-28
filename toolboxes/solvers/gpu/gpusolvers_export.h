/*
 * solvers_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GPUSOLVERS_EXPORT_H_
#define GPUSOLVERS_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPUSOLVERS__) || defined (gpusolvers_EXPORTS)
#define EXPORTGPUSOLVERS __declspec(dllexport)
#else
#define EXPORTGPUSOLVERS __declspec(dllimport)
#endif
#else
#define EXPORTGPUSOLVERS
#endif

#endif /* GPUSOLVERS_EXPORT_H_ */
