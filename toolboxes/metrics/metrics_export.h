#ifndef METRICS_EXPORT_H_
#define METRICS_EXPORT_H_

#if defined (WIN32)
   #if defined (__BUILD_GADGETRON_METRICS__) || defined (gadgetron_toolbox_METRICS__EXPORTS)
      #define EXPORTGADGETRONMETRICS __declspec(dllexport)
   #else
      #define EXPORTGADGETRONMETRICS __declspec(dllimport)
   #endif
#else
   #define EXPORTGADGETRONMETRICS
#endif

#endif /* METRICS_EXPORT_H_ */
