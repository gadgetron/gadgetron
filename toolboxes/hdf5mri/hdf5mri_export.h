/*
 * hdf5utils_export.h
 *
 *  Created on: Jan 20, 2012
 *      Author: Michael S. Hansen
 */

#ifndef HDF5MRI_EXPORT_H_
#define HDF5MRI_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_HDF5MRI__) || defined (hdf5mri_EXPORTS)
#define EXPORTHDF5MRI __declspec(dllexport)
#else
#define EXPORTHDF5MRI __declspec(dllimport)
#endif
#else
#define EXPORTHDF5MRI
#endif



#endif /* HDF5MRI_EXPORT_H_ */
