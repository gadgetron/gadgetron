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
	
	find_library(ZFP_DEBUG_LIBRARY NAMES ZFPd zfpd PATHS ${CMAKE_SOURCE_DIR}/../ZFP_wrappers/lib/ /usr/lib /usr/local/lib $ENV{ZFP_ROOT}/lib $ENV{ZFP_ROOT} DOC "ZFP library file (debug version)") 
endif () 
 
 
if (ZFP_INCLUDE_DIR AND ZFP_LIBRARY) 
	set(ZFP_FOUND TRUE) 
else () 
	set(ZFP_FOUND FALSE) 
endif () 
 
if (ZFP_DEBUG_LIBRARY) 
	set(ZFP_DEBUG_FOUND TRUE) 
else ()
  set(ZFP_DEBUG_LIBRARY ${ZFP_LIBRARY})
endif () 
 
if (ZFP_FOUND) 
	if (NOT ZFP_FIND_QUIETLY) 
		message(STATUS "Found ZFP library: ${ZFP_LIBRARY}") 
		message(STATUS "Found ZFP include: ${ZFP_INCLUDE_DIR}") 
	endif () 
else () 
	if (ZFP_FIND_REQUIRED) 
		message(FATAL_ERROR "Could not find ZFP") 
	endif () 
endif () 

# TSS: backwards compatibility
set(ZFP_LIBRARIES ${ZFP_LIBRARY}) 
