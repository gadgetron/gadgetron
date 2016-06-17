#ifndef FATWATER_EXPORT_H_
#define FATWATER_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_FATWATER__) || defined (fatwaterEXPORTS)
        #define EXPORTFATWATER __declspec(dllexport)
    #else
        #define EXPORTFATWATER __declspec(dllimport)
    #endif
#else
    #define EXPORTFATWATER
#endif

#endif // FATWATER_EXPORT_H_
