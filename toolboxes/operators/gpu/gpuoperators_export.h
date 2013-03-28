#ifndef GPUOPERATORS_EXPORT_H_
#define GPUOPERATORS_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPUOPERATORS__) || defined (gpusolvers_EXPORTS)
#define EXPORTGPUOPERATORS __declspec(dllexport)
#else
#define EXPORTGPUOPERATORS __declspec(dllimport)
#endif
#else
#define EXPORTGPUOPERATORS
#endif

#endif /* GPUOPERATORS_EXPORT_H_ */
