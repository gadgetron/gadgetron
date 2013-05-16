#ifndef GADGETRONSBSENSE_EXPORT_H_
#define GADGETRONSBSENSE_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_SBSENSE__) || defined (gadgetroncgsense_EXPORTS)
#define EXPORTGADGETSSBSENSE __declspec(dllexport)
#else
#define EXPORTGADGETSSBSENSE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSSBSENSE
#endif


#endif /* GADGETRONSBSENSE_EXPORT_H_ */
