#ifndef GADGETRON_HYPER_EXPORT_H_
#define GADGETRON_HYPER_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_HYPER__)
#define EXPORTGADGETSHYPER __declspec(dllexport)
#else
#define EXPORTGADGETSHYPER __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSHYPER
#endif

#endif /* GADGETRON_HYPER_EXPORT_H_ */
