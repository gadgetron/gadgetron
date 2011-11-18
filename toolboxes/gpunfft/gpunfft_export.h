/*
 * gpunfft_export.h
 *
 *  Created on: Nov 18, 2011
 *      Author: Michael S. Hansen
 */

#ifndef GPUNFFT_EXPORT_H_
#define GPUNFFT_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_GPUNFFT__) || define (gpunfft_EXPORTS)
#define EXPORTGPUNFFT __declspec(dllexport)
#else
#define EXPORTGPUNFFT __declspec(dllimport)
#endif
#else
#define EXPORTGPUNFFT
#endif


#endif /* GPUNFFT_EXPORT_H_ */
