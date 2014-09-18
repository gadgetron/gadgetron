#ifndef GADGETRON_GPUPMRI_EXPORT_H_
#define GADGETRON_GPUPMRI_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GADGET_GPUPMRI__)
#define EXPORTGADGETS_GPUPMRI __declspec(dllexport)
#else
#define EXPORTGADGETS_GPUPMRI __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_GPUPMRI
#endif

#endif /* GADGETRON_GPUPMRI_EXPORT_H_ */
