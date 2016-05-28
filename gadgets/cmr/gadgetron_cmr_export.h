#ifndef GADGETRON_CMR_EXPORT_H_
#define GADGETRON_CMR_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_CMR__)
        #define EXPORTGADGETSCMR __declspec(dllexport)
    #else
        #define EXPORTGADGETSCMR __declspec(dllimport)
    #endif
#else
    #define EXPORTGADGETSCMR
#endif

#endif /* GADGETRON_CMR_EXPORT_H_ */
