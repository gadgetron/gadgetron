# - Find the python 3 libraries
#   This module defines
#   PYTHONLIBS_FOUND           - have the Python libs been found
#   PYTHON_EXECUTABLE          - path to python executable
#   PYTHON_LIBRARIES           - path to the python library
#   PYTHON_INCLUDE_PATH        - path to where Python.h is found (deprecated)
#   PYTHON_INCLUDE_DIRS        - path to where Python.h is found
#   PYTHON_DEBUG_LIBRARIES     - path to the debug library (deprecated)
#   PYTHONLIBS_VERSION_STRING  - version of the Python libs found (since CMake 2.8.8)

message("Find python 3")

if (${CMAKE_VERSION}  VERSION_LESS "3.12.0")
    find_path(PYTHON3_PATH bin/python3 HINTS ENV{PYTHON3_PATH} PATHS /usr /usr/ /usr/local)
#


    find_path(PYTHON_INCLUDE_DIR3 NAMES python3.9/patchlevel.h python3.9m/patchlevel.h python3.7/patchlevel.h python3.7m/patchlevel.h python3.6/patchlevel.h python3.6m/patchlevel.h python3.5/patchlevel.h python3.5m/patchlevel.h  PATHS /usr/include /usr/local/include )
    if (EXISTS ${PYTHON_INCLUDE_DIR3}/python3.9m/patchlevel.h)
        set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIR3}/python3.9m)
    elseif    
    if (EXISTS ${PYTHON_INCLUDE_DIR3}/python3.7m/patchlevel.h)
        set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIR3}/python3.7m)
    elseif (EXISTS ${PYTHON_INCLUDE_DIR3}/python3.6m/patchlevel.h)
        set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIR3}/python3.6m)
    elseif (EXISTS ${PYTHON_INCLUDE_DIR3}/python3.5m/patchlevel.h)
        set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIR3}/python3.5m)
    endif()
            # find the python version
    if(EXISTS "${PYTHON_INCLUDE_DIRS}/patchlevel.h")
        file(STRINGS "${PYTHON_INCLUDE_DIRS}/patchlevel.h" python_version_str
            REGEX "^#define[ \t]+PY_VERSION[ \t]+\"[^\"]+\"")
        string(REGEX REPLACE "^#define[ \t]+PY_VERSION[ \t]+\"([^\"]+)\".*" "\\1"
                            PYTHONLIBS_VERSION_STRING "${python_version_str}")
        unset(python_version_str)
        message("Found python ${PYTHONLIBS_VERSION_STRING}")
    endif()

    string(REGEX MATCH "[0-9].[0-9]" PYTHON_MAJOR_VERSION ${PYTHONLIBS_VERSION_STRING})
    find_library(PYTHON_LIBRARIES libpython${PYTHON_MAJOR_VERSION}m.so)
    set(PYTHON_LIBRARY ${PYTHON_LIBRARIES})
    UNSET(PYTHON_EXECUTABLE CACHE)
    find_file(PYTHON_EXECUTABLE python3 PATHS /usr/bin /bin /usr/local/bin)
    set(PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIRS} )
    set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_DIRS})
    set(PYTHONLIBS_FOUND 1)

    else()
	if (Python3_compat_FIND_REQUIRED)
		find_package(Python3 COMPONENTS Interpreter REQUIRED)
	else()
		find_package(Python3 COMPONENTS Interpreter REQUIRED)
	endif()
	set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
	
endif ()

