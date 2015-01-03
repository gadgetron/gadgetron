#ifndef URLENCODE_H
#define URLENCODE_H

#include "log.h"

namespace Gadgetron {

/**
   Simple utility function for removing spaces and backslashes in URLs
   This function is used in various places to ensure proper encoding of schemalocation URIs
   
*/
inline std::string url_encode(const std::string& in) {
	char* tmp = new char[in.size()*4]; //Leave us plenty of space
	if (!tmp) {
		GDEBUG_STREAM("Failed to allocate temporary space for string in url_encode" << std::endl);
		return in;
	}

	char* optr = tmp;
	char* iptr = (char*)in.c_str();

	unsigned int counter = 0;
	while (counter < in.size()) {
		if (*iptr == ' ') {
			*optr++ = '%';
			*optr++ = '2';
			*optr++ = '0';
		} else if (*iptr == '\\') {
			*optr++ = '/';
		} else {
			*optr++ = *iptr;
		}
		iptr++;
		counter++;
	}
	*optr = '\0';

	std::string ret(tmp);

	delete [] tmp;

	return ret;
}
}

#endif //URLENCODE_H
