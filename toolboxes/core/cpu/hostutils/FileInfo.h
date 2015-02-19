#ifndef FILEINFO_H_
#define FILEINFO_H_

#include <string>
#include <fstream>

namespace Gadgetron {

/**
 *  Simple wrapper class for getting file info (file exists, file length, etc) before accessing file
 *
 */
class FileInfo
{
public:

	/**
	 *   Constructor. After calling the constructor the file_exists_ flag will be set in the class
	 */
	FileInfo(std::string filename)
	{
		filename_ = filename;
		std::ifstream ifile(filename_.c_str());
        file_exists_ = ifile.good();
	}

	virtual ~FileInfo() {}

	/**
	 *  Does the file exist (can be opened)
	 */
	bool exists() {
		return file_exists_;
	}

	size_t length() {
		size_t length = 0;
		if (file_exists_) {
			std::ifstream ifile(filename_.c_str());
			ifile.seekg(0,std::ios::end);
			length = ifile.tellg();
		} else {
			return -1;
		}
		return length;
	}

protected:
	bool file_exists_;
	std::string filename_;
};
}

#endif /* FILEINFO_H_ */
