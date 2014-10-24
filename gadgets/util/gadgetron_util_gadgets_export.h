#ifndef GADGETRON_UTIL_GADGETS_EXPORT_H_
#define GADGETRON_UTIL_GADGETS_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_UTIL_GADGETS__)
#define EXPORTUTILGADGETS __declspec(dllexport)
#else
#define EXPORTUTILGADGETS __declspec(dllimport)
#endif
#else
#define EXPORTUTILGADGETS
#endif

#endif // GADGETRON_UTIL_GADGETS_EXPORT_H_
