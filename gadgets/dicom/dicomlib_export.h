#ifndef DICOMLIB_EXPORT_H_
#define DICOMLIB_EXPORT_H_


#if defined (WIN32)
#if defined (gadgetrondicomlib_EXPORTS)
#define EXPORTGADGETSDICOMLIB __declspec(dllexport)
#else
#define EXPORTGADGETSDICOMLIB __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSDICOMLIB
#endif

#endif /* DICOMLIB_EXPORT_H_ */
