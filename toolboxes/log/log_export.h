#ifndef LOG_EXPORT_H_
#define LOG__EXPORT_H_

#if defined (WIN32)
    #ifdef BUILD_TOOLBOX_STATIC
        #define EXPORTGADGETRONLOG
    #else
        #if defined (__BUILD_GADGETRON_LOG___) || defined (gadgetron_toolbox_log__EXPORTS)
            #define EXPORTGADGETRONKLOG __declspec(dllexport)
        #else
            #define EXPORTGADGETRONLOG __declspec(dllimport)
        #endif
    #endif
#else
#define EXPORTGADGETRONLOG
#endif

#endif /* LOG__EXPORT_H_ */
