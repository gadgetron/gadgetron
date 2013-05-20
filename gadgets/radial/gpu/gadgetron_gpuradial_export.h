#ifndef GADGETRON_GPURADIAL_EXPORT_H_
#define GADGETRON_GPURADIAL_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPURADIAL__)
#define EXPORTGADGETS_GPURADIAL __declspec(dllexport)
#else
#define EXPORTGADGETS_GPURADIAL __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_GPURADIAL
#endif

#endif /* GADGETRON_GPURADIAL_EXPORT_H_ */
