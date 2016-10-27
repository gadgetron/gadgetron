#
# Find the Gadgetron Installation
#

# This module defines
# GADGETRON_INCLUDE_DIR, where to finds Gadget.h
# GADGETRON_HOME, Gadgetron Root Dir
# GADGETRON_LIB_DIR, This is where all the installed gadgetron libraries live
# GADGETRON_FOUND, if false, you cannot build anything that requires Gadgetron
# GADGETRON_VERSION_STRING, if Gadgetron is found, this verison string will be set

# Keep a list of variable names that we need to pass on to
# find_package_handle_standard_args().
set(_check_list)

# Search for the header file.
find_path(GADGETRON_HOME include/gadgetron/Gadget.h
    HINTS $ENV{GADGETRON_HOME} /usr/local/gadgetron /usr/gadgetron)
mark_as_advanced(GADGETRON_HOME)
list(APPEND _check_list GADGETRON_HOME)

set(GADGETRON_INCLUDE_DIR ${GADGETRON_HOME}/include/gadgetron)
mark_as_advanced(GADGETRON_INCLUDE_DIR)
list(APPEND _check_list GADGETRON_INCLUDE_DIR)

set(GADGETRON_LIB_DIR ${GADGETRON_HOME}/lib)
mark_as_advanced(GADGETRON_LIB_DIR)
list(APPEND _check_list GADGETRON_LIB_DIR)

# Handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gadgetron DEFAULT_MSG ${_check_list})

# If Cuda is detected on the system some header files will be needed
# -- whether Cuda is used or not --

find_package(CUDA)
if (CUDA_FOUND)
  include_directories( ${CUDA_INCLUDE_DIRS} )
endif ()

# ------------------------------------------------------------------------
#  Extract version information from "incluce/gadgetron/gadgetron_config.h"
# ------------------------------------------------------------------------
if(GADGETRON_FOUND)
    set(GADGETRON_VERSION_MAJOR 0)
    set(GADGETRON_VERSION_MINOR 0)
    set(GADGETRON_VERSION_PATCH 0) 

    if(EXISTS "${GADGETRON_INCLUDE_DIR}/gadgetron_config.h")

        # Read and parse header file for version number
        file(STRINGS "${GADGETRON_INCLUDE_DIR}/gadgetron_config.h" _gadgetron_HEADER_CONTENTS REGEX "#define GADGETRON_VERSION_MAJOR ")
        string(REGEX REPLACE ".*#define GADGETRON_VERSION_MAJOR ([0-9]+).*" "\\1" GADGETRON_VERSION_MAJOR "${_gadgetron_HEADER_CONTENTS}")
        unset(_gadgetron_HEADER_CONTENTS)
        
        file(STRINGS "${GADGETRON_INCLUDE_DIR}/gadgetron_config.h" _gadgetron_HEADER_CONTENTS REGEX "#define GADGETRON_VERSION_MINOR ")
        string(REGEX REPLACE ".*#define GADGETRON_VERSION_MINOR ([0-9]+).*" "\\1" GADGETRON_VERSION_MINOR "${_gadgetron_HEADER_CONTENTS}")
        unset(_gadgetron_HEADER_CONTENTS)
        
        file(STRINGS "${GADGETRON_INCLUDE_DIR}/gadgetron_config.h" _gadgetron_HEADER_CONTENTS REGEX "#define GADGETRON_VERSION_PATCH ")
        string(REGEX REPLACE ".*#define GADGETRON_VERSION_PATCH ([0-9]+).*" "\\1" GADGETRON_VERSION_PATCH "${_gadgetron_HEADER_CONTENTS}")
        unset(_gadgetron_HEADER_CONTENTS)
        
    endif (EXISTS "${GADGETRON_INCLUDE_DIR}/gadgetron_config.h")

    set(GADGETRON_VERSION_STRING ${GADGETRON_VERSION_MAJOR}.${GADGETRON_VERSION_MINOR}.${GADGETRON_VERSION_PATCH})
    set(GADGETRON_SOVERSION ${GADGETRON_VERSION_MAJOR}.${GADGETRON_VERSION_MINOR})
    message("find Gadgetron version : ${GADGETRON_VERSION_STRING}")
endif()