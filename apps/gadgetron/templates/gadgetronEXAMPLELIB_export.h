/*
 * gadgetronEXAMPLELIB_export.h
 *
 */

#ifndef GADGETRONEXAMPLELIB_EXPORT_H_
#define GADGETRONEXAMPLELIB_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_EXAMPLELIB__) || defined (gadgetronEXAMPLELIB_EXPORTS)
#define EXPORTGADGETSEXAMPLELIB __declspec(dllexport)
#else
#define EXPORTGADGETSEXAMPLELIB __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSEXAMPLELIB
#endif


#endif /* GADGETRONEXAMPLELIB_EXPORT_H_ */
