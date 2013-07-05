# - this module looks for Matlab
# Defines:
#  MATLAB_INCLUDE_DIR:  include path for mex.h, engine.h
#  MATLAB_LIBRARIES:    required libraries: libmex, etc
#  MATLAB_JARS:         optional java jars: jmi.jar, util.jar, etc
#  MATLAB_MEX_LIBRARY:  path to libmex.lib
#  MATLAB_MX_LIBRARY:   path to libmx.lib
#  MATLAB_ENG_LIBRARY:  path to libeng.lib

#=============================================================================
# Copyright 2005-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)
#
#  Modified April 2013 - Joseph Naegele
#   - Updated to work on OS X 10._
#   - Added Matlab's Java Jars as MATLAB_JARS

set(MATLAB_FOUND 0)

if(WIN32)
    set(MATLAB_ROOT "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.0;MATLABROOT]")
    if(${CMAKE_GENERATOR} MATCHES "Visual Studio 6")
        set(MATLAB_LIB_ROOT "${MATLAB_ROOT}/extern/lib/win32/microsoft/msvc60")
    else()
        if(${CMAKE_GENERATOR} MATCHES "Visual Studio 7")
            # Assume people are generally using 7.1,
            # if using 7.0 need to link to: ../extern/lib/win32/microsoft/msvc70
            set(MATLAB_LIB_ROOT "${MATLAB_ROOT}/extern/lib/win32/microsoft/msvc71")
        else()
            if(${CMAKE_GENERATOR} MATCHES "Borland")
                # Same here, there are also: bcc50 and bcc51 directories
                set(MATLAB_LIB_ROOT "${MATLAB_ROOT}/extern/lib/win32/microsoft/bcc54")
            else()
                if(MATLAB_FIND_REQUIRED)
                    message(FATAL_ERROR "Generator not compatible: ${CMAKE_GENERATOR}")
                endif()
            endif()
        endif()
    endif()
    find_path(
        MATLAB_INCLUDE_DIR
        "mex.h"
        HINTS ${MATLAB_ROOT}
        PATH_SUFFIXES extern/include
    )
else(WIN32)
    find_path(
        MATLAB_ROOT extern/include/mex.h
        HINTS $ENV{MATLAB_HOME} $ENV{MATLAB_ROOT}
        PATHS /usr /usr/local /opt
        PATH_SUFFIXES MATLAB
    )
    find_path(
        MATLAB_INCLUDE_DIR
        "mex.h"
        HINTS ${MATLAB_ROOT}
        PATH_SUFFIXES extern/include
    )
endif(WIN32)

# find each library
# give it it's own cmake variable
# add it to the list of libraries
foreach(lib mat ut mex mx eng)
    string(TOUPPER ${lib} LIB)
    find_library(
        MATLAB_${LIB}_LIBRARY
        ${lib}
        HINTS ${MATLAB_ROOT} ${MATLAB_LIB_ROOT}
        PATH_SUFFIXES lib bin bin/maci64 bin/glnxa64 bin/glnax86
    )
    if(MATLAB_${LIB}_LIBRARY)
        list(APPEND MATLAB_LIBRARIES "${MATLAB_${LIB}_LIBRARY}")
    endif()
    mark_as_advanced(MATLAB_${LIB}_LIBRARY)
endforeach()

foreach(jar jmi util)
    string(TOUPPER ${jar} LIB)
    find_file(
        MATLAB_${LIB}_JAR
        "${jar}.jar"
        HINTS ${MATLAB_ROOT}
        PATH_SUFFIXES java jar java/jar
    )
    if(MATLAB_${LIB}_JAR)
        list(APPEND MATLAB_JARS "${MATLAB_${LIB}_JAR}")
    endif()
    mark_as_advanced(MATLAB_${LIB}_JAR)
endforeach()

if(MATLAB_INCLUDE_DIR AND MATLAB_LIBRARIES)
    set(MATLAB_FOUND 1)
endif()

include("FindPackageHandleStandardArgs")
FIND_PACKAGE_HANDLE_STANDARD_ARGS("Matlab" DEFAULT_MSG MATLAB_ROOT MATLAB_INCLUDE_DIR MATLAB_LIBRARIES)

mark_as_advanced(
    MATLAB_JARS
    MATLAB_LIBRARIES
    MATLAB_INCLUDE_DIR
    MATLAB_FOUND
    MATLAB_ROOT
)

