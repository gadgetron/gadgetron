############################################################################
# FindBART
#
# Find BART includes and libraries
# 
# This module defines the following variables:
#
# BART_ROOT            -where to find bart root 
# BART_INCLUDE_DIR     -where to find the main.h and bart_embed_api.h
# BART_LIBARIES_DIR    -List of libraries when using bart
# BART_DIR             -where to find bart source code 
# BART_FOUND           -Set to True if bart is found
#
# BART_VERSION_STRING  -The version of bart found (major.minor.patch)
# BART_VERSION_MAJOR   -The major version of bart
# BART_VERSION_MINOR   -The minor version of bart
# BART_VERSION_PATCH   -The patch version of bart
############################################################################

# Set default directories
set(BART_ROOT /usr/local)
set(BART_INCLUDE_DIR ${BART_ROOT}/include/bart)
set(BART_LIBRARIES_DIR ${BART_ROOT}/lib)
set(BART_DIR /opt/code/bart)

set(BART_VERSION_STRING "")
set(BART_VERSION_MAJOR 0)
set(BART_VERSION_MINOR 0)
set(BART_VERSION_PATCH 0) 

find_path(BART_INCLUDE_DIR 
		NAMES
		main.h
		bart_embed_api.h
		PATHS
		${BART_INCLUDE_DIR})

find_library(BART_LIBRARY
		NAMES
		bartmain
		bartsupport
		PATHS
		${BART_LIBRARIES_DIR})

if(BART_INCLUDE_DIR AND BART_LIBRARY) 
	set(BART_FOUND TRUE) 
else() 
	set(BART_FOUND FALSE) 
endif() 

# Get BART version
if(BART_FOUND AND (EXISTS ${BART_DIR}/.git))
  find_package(Git QUIET)
  if(Git_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --match "v*" --dirty
      OUTPUT_VARIABLE full_version
      RESULT_VARIABLE ret
      WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})

    if(ret EQUAL 0)
      string(STRIP ${full_version} full_version)
      string(REGEX MATCHALL "[0-9]+" version_list ${full_version})
      list(GET version_list 0 BART_VERSION_MAJOR)
      list(GET version_list 1 BART_VERSION_MINOR)
      list(GET version_list 2 BART_VERSION_PATCH)
      set(BART_VERSION_STRING ${BART_VERSION_MAJOR}.${BART_VERSION_MINOR}.${BART_VERSION_PATCH})     
      message("find BART version : ${BART_VERSION_STRING}")
    else()
      message("find BART")
    endif()

 endif()

endif()

