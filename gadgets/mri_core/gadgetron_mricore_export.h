#ifndef GADGETRON_MRICORE_EXPORT_H_
#define GADGETRON_MRICORE_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_MRICORE__)
#define EXPORTGADGETSMRICORE __declspec(dllexport)
#else
#define EXPORTGADGETSMRICORE __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSMRICORE
#endif

#endif /* GADGETRON_MRICORE_EXPORT_H_ */
