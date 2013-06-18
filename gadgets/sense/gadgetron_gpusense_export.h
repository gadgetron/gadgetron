#ifndef GADGETRON_GPUSENSE_EXPORT_H_
#define GADGETRON_GPUSENSE_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPUSENSE__)
#define EXPORTGADGETS_GPUSENSE __declspec(dllexport)
#else
#define EXPORTGADGETS_GPUSENSE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_GPUSENSE
#endif

#endif /* GADGETRON_GPUSENSE_EXPORT_H_ */
