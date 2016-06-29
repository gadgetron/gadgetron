#
# Find the ZFP includes and library
#

# This module defines
# ZFP_INCLUDE_DIR, where to find zfp.h
# ZFP_LIBRARIES, the libraries to link against
# ZFP_FOUND, if false, you cannot build anything that requires ZFP

######################################################################## 
find_path(ZFP_INCLUDE_DIR zfp/zfp.h /usr/include /usr/local/include $ENV{ZFP_ROOT} $ENV{ZFP_ROOT}/inc DOC "directory containing zfp/zfp.h for ZFP library")  
find_library(ZFP_LIBRARY NAMES ZFP zfp PATHS /usr/lib /usr/local/lib $ENV{ZFP_ROOT}/lib $ENV{ZFP_ROOT} DOC "ZFP library file") 
 
if (WIN32 AND NOT CYGWIN) 
	set(CMAKE_DEBUG_POSTFIX "d") 
	find_library(ZFP_DEBUG_LIBRARY NAMES ZFP${CMAKE_DEBUG_POSTFIX} zfp${CMAKE_DEBUG_POSTFIX} PATHS ${CMAKE_SOURCE_DIR}/../ZFP_wrappers/lib/ /usr/lib /usr/local/lib $ENV{ZFP_ROOT}/lib $ENV{ZFP_ROOT} DOC "ZFP library file (debug version)") 
endif (WIN32 AND NOT CYGWIN) 
 
 
if (ZFP_INCLUDE_DIR AND ZFP_LIBRARY) 
	set(ZFP_FOUND TRUE) 
else (ZFP_INCLUDE_DIR AND ZFP_LIBRARY) 
	set(ZFP_FOUND FALSE) 
endif (ZFP_INCLUDE_DIR AND ZFP_LIBRARY) 
 
if (ZFP_DEBUG_LIBRARY) 
	set(ZFP_DEBUG_FOUND TRUE) 
else (ZFP_DEBUG_LIBRARY)
  set(ZFP_DEBUG_LIBRARY ${ZFP_LIBRARY})
endif (ZFP_DEBUG_LIBRARY) 
 
if (ZFP_FOUND) 
	if (NOT ZFP_FIND_QUIETLY) 
		message(STATUS "Found ZFP library: ${ZFP_LIBRARY}") 
		message(STATUS "Found ZFP include: ${ZFP_INCLUDE_DIR}") 
	endif (NOT ZFP_FIND_QUIETLY) 
else (ZFP_FOUND) 
	if (ZFP_FIND_REQUIRED) 
		message(FATAL_ERROR "Could not find ZFP") 
	endif (ZFP_FIND_REQUIRED) 
endif (ZFP_FOUND) 

# TSS: backwards compatibility
set(ZFP_LIBRARIES ${ZFP_LIBRARY}) 
