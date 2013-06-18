#ifndef GADGETRON_RADIAL_EXPORT_H_
#define GADGETRON_RADIAL_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_RADIAL__)
#define EXPORTGADGETS_RADIAL __declspec(dllexport)
#else
#define EXPORTGADGETS_RADIAL __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_RADIAL
#endif

#endif /* GADGETRON_GPURADIAL_EXPORT_H_ */
