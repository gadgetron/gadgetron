/** \file cpunfft_export.h
    \brief Required definitions for Windows, importing/exporting dll symbols 
*/

#ifndef CPUNFFT_EXPORT_H_
#define CPUNFFT_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_CPUNFFT__) || defined (cpunfft_EXPORTS)
#define EXPORTCPUNFFT __declspec(dllexport)
#else
#define EXPORTCPUNFFT __declspec(dllimport)
#endif
#else
#define EXPORTCPUNFFT
#endif


#endif /* CPUNFFT_EXPORT_H_ */
