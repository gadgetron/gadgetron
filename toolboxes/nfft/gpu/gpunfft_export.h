/** \file gpunfft_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef GPUNFFT_EXPORT_H_
#define GPUNFFT_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPUNFFT__) || defined (gpunfft_EXPORTS)
#define EXPORTGPUNFFT __declspec(dllexport)
#else
#define EXPORTGPUNFFT __declspec(dllimport)
#endif
#else
#define EXPORTGPUNFFT
#endif


#endif /* GPUNFFT_EXPORT_H_ */
