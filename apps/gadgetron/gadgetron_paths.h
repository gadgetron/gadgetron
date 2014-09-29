#ifndef GADGETRON_PATHS_H
#define GADGETRON_PATHS_H

#include <limits.h>
#include <unistd.h>
#include <string>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#include <Shlwapi.h>
#pragma comment(lib, "shlwapi.lib")
#endif // _WIN32

#ifdef __APPLE__
#include <mach-o/dyld.h>/* _NSGetExecutablePath */
#endif

#define MAX_GADGETRON_HOME_LENGTH 1024

namespace Gadgetron
{
  inline std::string get_gadgetron_home()
  {
#if defined  __APPLE__
    char path[MAX_GADGETRON_HOME_LENGTH];
    uint32_t size = sizeof(path);
    if (_NSGetExecutablePath(path, &size) == 0) {
      std::string s1(path);
      return s1.substr(0, s1.find_last_of("\\/")) + std::string("/../");
    } else {
      std::cout << "Unable to determine GADGETRON_HOME" << std::endl;
      return std::string("");
    }
#elif defined _WIN32 || _WIN64
    // Full path to the executable (including the executable file)
    char fullPath[MAX_GADGETRON_HOME_LENGTH];	
    // Full path to the executable (without executable file)
    char *rightPath;
    // Will contain exe path
    HMODULE hModule = GetModuleHandle(NULL);
    if (hModule != NULL)
      {
	// When passing NULL to GetModuleHandle, it returns handle of exe itself
	GetModuleFileName(hModule, fullPath, (sizeof(fullPath))); 
	rightPath = fullPath;
	PathRemoveFileSpec(rightPath);
	std::string s1(rightPath);
	return s1 + std::string("\\..\\");
      }
    else
      {
        std::cout << "The path to the executable is NULL" << std::endl;
      }
#else //Probably some NIX where readlink should work
    char buff[MAX_GADGETRON_HOME_LENGTH];
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1) {
      buff[len] = '\0';
      std::string s1(buff);
      return s1.substr(0, s1.find_last_of("\\/")) + std::string("/../");
    } else {
      std::cout << "Unable to determine GADGETRON_HOME" << std::endl;
      return std::string("");
    }
#endif
  }
}

#endif //GADGETRON_PATHS_H
