/** \file cpucore_math_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
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
