/*
 * hdf5utils_export.h
 *
 *  Created on: Jan 20, 2012
 *      Author: Michael S. Hansen
 */

#ifndef HDF5UTILS_EXPORT_H_
#define HDF5UTILS_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_HDF5UTILS__) || defined (hdf5utils_EXPORTS)
#define EXPORTHDF5UTILS __declspec(dllexport)
#else
#define EXPORTHDF5UTILS __declspec(dllimport)
#endif
#else
#define EXPORTHDF5UTILS
#endif



#endif /* HDF5UTILS_EXPORT_H_ */
