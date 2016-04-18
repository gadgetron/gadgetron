#ifndef GADGETRON_FATWATER_EXPORT_H_
#define GADGETRON_FATWATER_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETS_FATWATER__) || defined (gadgetron_fatwaterEXPORTS)
        #define EXPORTGADGETFATWATER __declspec(dllexport)
    #else
        #define EXPORTGADGETFATWATER __declspec(dllimport)
    #endif
#else
    #define EXPORTGADGETFATWATER
#endif

#endif // GADGETRON_FATWATER_EXPORT_H_
