/** \file nfft_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef NFFT_EXPORT_H_
#define NFFT_EXPORT_H_

#if defined (WIN32)
    #if defined (__BUILD_GADGETRON_NFFT__) || defined (fft_EXPORTS)
        #define EXPORTNFFT __declspec(dllexport)
    #else
        #define EXPORTNFFT __declspec(dllimport)
    #endif
#else
    #define EXPORTNFFT
#endif


#endif /* NFFT_EXPORT_H_ */
