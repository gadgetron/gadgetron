#ifndef _GPUREG_EXPORT_H_
#define _GPUREG_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPU__) || defined (gpureg_EXPORTS)
#define EXPORTGPUREG __declspec(dllexport)
#else
#define EXPORTGPUREG __declspec(dllimport)
#endif
#else
#define EXPORTGPUREG
#endif

#endif /* _GPUREG_EXPORT_H_ */
