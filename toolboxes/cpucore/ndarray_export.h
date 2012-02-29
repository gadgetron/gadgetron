/*
 * ndarray_export.h
 *
 *  Created on: Nov 29, 2011
 *      Author: hansenms
 */

#ifndef NDARRAY_EXPORT_H_
#define NDARRAY_EXPORT_H_


#if defined (WIN32)
#if defined (__BUILD_GADGETRON_NDARRAY__) || defined (hondarray_EXPORTS)
#define EXPORTNDARRAY __declspec(dllexport)
#else
#define EXPORTNDARRAY __declspec(dllimport)
#endif
#else
#define EXPORTNDARRAY
#endif



#endif /* NDARRAY_EXPORT_H_ */
