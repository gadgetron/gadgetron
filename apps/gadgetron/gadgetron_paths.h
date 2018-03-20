#ifndef GADGETRON_PATHS_H
#define GADGETRON_PATHS_H

#include <limits.h>
#include <string>
#include "log.h"

#ifdef _WIN32
#include <windows.h>
#include <Shlwapi.h>
#pragma comment(lib, "shlwapi.lib")
#else
#include <unistd.h>
#include <cstdlib>
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
    char path[PATH_MAX];
    uint32_t size = sizeof(path);
    char resolved[PATH_MAX];
    if ((_NSGetExecutablePath(path, &size) == 0) && (realpath(path, resolved) != NULL)) {
      std::string s1(resolved);
      return s1.substr(0, s1.find_last_of("\\/")) + std::string("/../");
    } else {
      GDEBUG_STREAM("Unable to determine GADGETRON_HOME" << std::endl);
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
	for(int i = 0; i < strlen(rightPath); i++)
	  if(rightPath[i] == '\\') rightPath[i] = '/';

	std::string s1(rightPath);
	return s1 + std::string("/../");
      }
    else
      {
        GDEBUG_STREAM("The path to the executable is NULL" << std::endl);
        return std::string("");
      }
#else //Probably some NIX where readlink should work
      std::string home = std::getenv("GADGETRON_HOME");

    if (home.size() == 0) {
      GDEBUG_STREAM("Unable to determine GADGETRON_HOME" << std::endl);
      return std::string("");
    } else {
        return home;
    }

#endif
  }
}

#endif //GADGETRON_PATHS_H
