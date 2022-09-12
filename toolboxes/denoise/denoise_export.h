#ifndef DENOISE_EXPORT_H_
#define DENOISE_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_DENOISE__) || defined (fatwaterEXPORTS)
        #define EXPORTDENOISE __declspec(dllexport)
    #else
        #define EXPORTDENOISE __declspec(dllimport)
    #endif
#else
    #define EXPORTDENOISE
#endif

#endif // DENOISE_EXPORT_H_
