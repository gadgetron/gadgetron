/*
 * cpucore_math_export.h
 *
 *  Created on: Nov 29, 2011
 *      Author: hansenms
 */

#ifndef CPUCORE_MATH_EXPORT_H_
#define CPUCORE_MATH_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_CPUCORE_MATH__) || defined (cpucore_math_EXPORTS)
#define EXPORTCPUCOREMATH __declspec(dllexport)
#else
#define EXPORTCPUCOREMATH __declspec(dllimport)
#endif
#else
#define EXPORTCPUCOREMATH
#endif

#endif /* CPUCORE_MATH_EXPORT_H_ */
