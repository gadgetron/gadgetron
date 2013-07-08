#ifndef GADGETRON_GRAPPA_EXPORT_H_
#define GADGETRON_GRAPPA_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GRAPPA__)
#define EXPORTGADGETSGRAPPA __declspec(dllexport)
#else
#define EXPORTGADGETSGRAPPA __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSGRAPPA
#endif

#endif /* GADGETRON_GRAPPA_EXPORT_H_ */
