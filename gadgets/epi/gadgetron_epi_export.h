#ifndef GADGETRON_EPI_EXPORT_H_
#define GADGETRON_EPI_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_EPI__)
#define EXPORTGADGETS_EPI __declspec(dllexport)
#else
#define EXPORTGADGETS_EPI __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_EPI
#endif

#endif /* GADGETRON_EPI_EXPORT_H_ */
