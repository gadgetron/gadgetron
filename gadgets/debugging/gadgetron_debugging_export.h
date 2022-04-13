#ifndef GADGETRON_DEBUGGING_EXPORT_H_
#define GADGETRON_DEBUGGING_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_DEBUGGING__)
#define EXPORTGADGETSDEBUGGING __declspec(dllexport)
#else
#define EXPORTGADGETSDEBUGGING __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSDEBUGGING
#endif

#endif /* GADGETRON_DEBUGGING_EXPORT_H_ */
