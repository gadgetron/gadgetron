#
# Find the ZFP includes and library
#

# This module defines
# ZFP_INCLUDE_DIR, where to find zfp.h
# ZFP_LIBRARIES, the libraries to link against
# ZFP_FOUND, if false, you cannot build anything that requires ZFP

######################################################################## 
FIND_PATH(ZFP_INCLUDE_DIR zfp/zfp.h /usr/include /usr/local/include $ENV{ZFP_ROOT} $ENV{ZFP_ROOT}/inc DOC "directory containing zfp/zfp.h for ZFP library")  
FIND_LIBRARY(ZFP_LIBRARY NAMES ZFP zfp PATHS /usr/lib /usr/local/lib $ENV{ZFP_ROOT}/lib $ENV{ZFP_ROOT} DOC "ZFP library file") 
 
IF (WIN32 AND NOT CYGWIN) 
	SET(CMAKE_DEBUG_POSTFIX "d") 
	FIND_LIBRARY(ZFP_DEBUG_LIBRARY NAMES ZFP${CMAKE_DEBUG_POSTFIX} zfp${CMAKE_DEBUG_POSTFIX} PATHS ${CMAKE_SOURCE_DIR}/../ZFP_wrappers/lib/ /usr/lib /usr/local/lib $ENV{ZFP_ROOT}/lib $ENV{ZFP_ROOT} DOC "ZFP library file (debug version)") 
ENDIF (WIN32 AND NOT CYGWIN) 
 
 
IF (ZFP_INCLUDE_DIR AND ZFP_LIBRARY) 
	SET(ZFP_FOUND TRUE) 
ELSE (ZFP_INCLUDE_DIR AND ZFP_LIBRARY) 
	SET(ZFP_FOUND FALSE) 
ENDIF (ZFP_INCLUDE_DIR AND ZFP_LIBRARY) 
 
IF (ZFP_DEBUG_LIBRARY) 
	SET(ZFP_DEBUG_FOUND TRUE) 
ELSE (ZFP_DEBUG_LIBRARY)
  SET(ZFP_DEBUG_LIBRARY ${ZFP_LIBRARY})
ENDIF (ZFP_DEBUG_LIBRARY) 
 
IF (ZFP_FOUND) 
	IF (NOT ZFP_FIND_QUIETLY) 
		MESSAGE(STATUS "Found ZFP library: ${ZFP_LIBRARY}") 
		MESSAGE(STATUS "Found ZFP include: ${ZFP_INCLUDE_DIR}") 
	ENDIF (NOT ZFP_FIND_QUIETLY) 
ELSE (ZFP_FOUND) 
	IF (ZFP_FIND_REQUIRED) 
		MESSAGE(FATAL_ERROR "Could not find ZFP") 
	ENDIF (ZFP_FIND_REQUIRED) 
ENDIF (ZFP_FOUND) 

# TSS: backwards compatibility
SET(ZFP_LIBRARIES ${ZFP_LIBRARY}) 
