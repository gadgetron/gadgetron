/*
 * hdf5_core.h
 *
 *  Created on: Jan 26, 2012
 *      Author: Michael S. Hansen
 */

#ifndef HDF5_CORE_H_
#define HDF5_CORE_H_

#include "hdf5utils_export.h"

#include <ace/Synch.h>
#include <ace/Mutex.h>

#include <FileInfo.h>
#include <H5Cpp.h>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

#ifndef H5_NO_NAMESPACE
	using namespace H5;
#endif


EXPORTHDF5UTILS boost::shared_ptr<H5File> OpenHDF5File(const char* filename);

EXPORTHDF5UTILS bool HDF5LinkExists(H5File* f, const char* name);

EXPORTHDF5UTILS int HDF5CreateGroupForDataset(H5File* f, const char* name);

EXPORTHDF5UTILS unsigned long HDF5GetLengthOfFirstDimension(const char* filename, const char* name);


class EXPORTHDF5UTILS HDF5Lock
{

public:
	static HDF5Lock* instance();

	void acquire();
	void release();

protected:
	HDF5Lock()
	: mutex_("HDF5ThreadMutex") { }

	virtual ~HDF5Lock() { }

	static HDF5Lock* instance_;

	ACE_Thread_Mutex mutex_;
};

class HDF5Exclusive
{
public:
	HDF5Exclusive() {
		HDF5Lock::instance()->acquire();
	}

	~HDF5Exclusive() {
		HDF5Lock::instance()->release();
	}

};

#endif /* HDF5_CORE_H_ */
