#ifndef LOG_EXPORT_H_
#define LOG_EXPORT_H_

#if defined (WIN32)
   #if defined (__BUILD_GADGETRON_LOG___) || defined (gadgetron_toolbox_log__EXPORTS)
      #define EXPORTGADGETRONKLOG __declspec(dllexport)
   #else
      #define EXPORTGADGETRONLOG __declspec(dllimport)
   #endif
#else
   #define EXPORTGADGETRONLOG
#endif

#endif /* LOG_EXPORT_H_ */
