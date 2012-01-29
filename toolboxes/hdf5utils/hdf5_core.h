/*
 * hdf5_core.h
 *
 *  Created on: Jan 26, 2012
 *      Author: Michael S. Hansen
 */

#ifndef HDF5_CORE_H_
#define HDF5_CORE_H_

#include "hdf5utils_export.h"

#include <FileInfo.h>
#include <H5Cpp.h>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

#ifndef H5_NO_NAMESPACE
	using namespace H5;
#endif


EXPORTHDF5UTILS  boost::shared_ptr<H5File> OpenHF5File(const char* filename);

EXPORTHDF5UTILS  bool HDF5LinkExists(H5File* f, const char* name);

EXPORTHDF5UTILS int HDF5CreateGroupForDataset(H5File* f, const char* name);

#endif /* HDF5_CORE_H_ */
