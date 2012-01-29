/*
 * hdf5_core.cpp
 *
 *  Created on: Jan 26, 2012
 *      Author: Michael S. Hansen
 */

#include "hdf5_core.h"


#ifndef H5_NO_NAMESPACE
	using namespace H5;
#endif

#include <vector>

boost::shared_ptr<H5File> OpenHF5File(const char* filename)
{
	boost::shared_ptr<H5File> ret;

	if (FileInfo(std::string(filename)).exists()) {
		ret = boost::shared_ptr<H5::H5File>(new H5File(filename, H5F_ACC_RDWR));
	} else {
		ret = boost::shared_ptr<H5::H5File>(new H5File(filename, H5F_ACC_TRUNC));
	}
	return ret;
}

bool HDF5LinkExists(H5File* f, const char* name)
{
	std::vector<std::string> name_elements;
	std::string splitstr("/");
	std::string namestr(name);
	boost::split(name_elements, namestr, boost::is_any_of(splitstr));
	std::string current_path("");
	for (unsigned int i = 0; i < name_elements.size(); i++) {
		if (name_elements[i].size() > 0) {
			current_path = current_path + std::string("/") + name_elements[i];
			if (!H5Lexists(f->getId(), current_path.c_str(), H5P_DEFAULT )) {
				return false;
			}
		}
	}
	return true;
}

int HDF5CreateGroupForDataset(H5File* f, const char* name)
{
	std::vector<std::string> name_elements;
	std::string splitstr("/");
	std::string namestr(name);
	boost::split(name_elements, namestr, boost::is_any_of(splitstr));
	std::string current_path("");
	for (unsigned int i = 0; i < name_elements.size()-1; i++) { //Skip the last one, we are just creating the group
		if (name_elements[i].size() > 0) {
			current_path = current_path + std::string("/") + name_elements[i];
			if (!H5Lexists(f->getId(), current_path.c_str(), H5P_DEFAULT )) {
				f->createGroup( current_path.c_str());
			}
		}
	}
	return 0;
}




