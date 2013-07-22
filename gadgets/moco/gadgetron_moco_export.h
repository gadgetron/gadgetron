#ifndef GADGETRON_MOCO_EXPORT_H_
#define GADGETRON_MOCO_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_MOCO__)
#define EXPORTGADGETS_MOCO __declspec(dllexport)
#else
#define EXPORTGADGETS_MOCO __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_MOCO
#endif

#endif /* GADGETRON_MOCO_EXPORT_H_ */
