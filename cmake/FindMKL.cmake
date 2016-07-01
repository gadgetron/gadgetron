# - Find the MKL libraries
# Modified from Armadillo's ARMA_FindMKL.cmake
# This module defines
#  MKL_INCLUDE_DIR, the directory for the MKL headers
#  MKL_LIB_DIR, the directory for the MKL library files
#  MKL_COMPILER_LIB_DIR, the directory for the MKL compiler library files
#  MKL_LIBRARIES, the libraries needed to use Intel's implementation of BLAS & LAPACK.
#  MKL_FOUND, If false, do not try to use MKL; if true, the macro definition USE_MKL is added.

# Set the include path
# TODO: what if MKL is not installed in /opt/intel/mkl?
# try to find at /opt/intel/mkl
# in windows, try to find MKL at C:/Program Files (x86)/Intel/Composer XE/mkl

if ( WIN32 )
  if(NOT DEFINED ENV{MKLROOT_PATH})
    set(MKLROOT_PATH "C:/Program Files (x86)/Intel/Composer XE" CACHE PATH "Where the MKL are stored")
  else(NOT DEFINED ENV{MKLROOT_PATH})
    message("MKLROOT_PATH is set as the corresponding environmental variable ... ")
    set(MKLROOT_PATH $ENV{MKLROOT_PATH} CACHE PATH "Where the MKL are stored")
  endif(NOT DEFINED ENV{MKLROOT_PATH}) 
else ( WIN32 )
    set(MKLROOT_PATH "/opt/intel" CACHE PATH "Where the MKL are stored")
endif ( WIN32 )

if (EXISTS ${MKLROOT_PATH}/mkl)
    SET(MKL_FOUND TRUE)
    message("MKL is found at ${MKLROOT_PATH}/mkl")
    IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set( USE_MKL_64BIT On )
        if ( ARMADILLO_FOUND )
            if ( ARMADILLO_BLAS_LONG_LONG )
                set( USE_MKL_64BIT_LIB On )
                ADD_DEFINITIONS(-DMKL_ILP64)
                message("MKL is linked against ILP64 interface ... ")
            endif ( ARMADILLO_BLAS_LONG_LONG )
        endif ( ARMADILLO_FOUND )
    ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set( USE_MKL_64BIT Off )
    ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8)
else (EXISTS ${MKLROOT_PATH}/mkl)
    SET(MKL_FOUND FALSE)
    message("MKL is NOT found ... ")
endif (EXISTS ${MKLROOT_PATH}/mkl)

if (MKL_FOUND)
    set(MKL_INCLUDE_DIR "${MKLROOT_PATH}/mkl/include")
    ADD_DEFINITIONS(-DUSE_MKL)
    if ( USE_MKL_64BIT )
        set(MKL_LIB_DIR "${MKLROOT_PATH}/mkl/lib/intel64")
        set(MKL_COMPILER_LIB_DIR "${MKLROOT_PATH}/compiler/lib/intel64")
        set(MKL_COMPILER_LIB_DIR ${MKL_COMPILER_LIB_DIR} "${MKLROOT_PATH}/lib/intel64")
        if (NOT WIN32)
            if ( USE_MKL_64BIT_LIB )
                    if (WIN32)
                        set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_ilp64_dll)
                    else (WIN32)
                        set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_ilp64)
                    endif (WIN32)
            else ( USE_MKL_64BIT_LIB )
                    if (WIN32)
                        set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_lp64_dll)
                    else (WIN32)
                        set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_lp64)
                    endif (WIN32)
            endif ( USE_MKL_64BIT_LIB )
        endif (NOT WIN32)
    else ( USE_MKL_64BIT )
        set(MKL_LIB_DIR "${MKLROOT_PATH}/mkl/lib/ia32")
        set(MKL_COMPILER_LIB_DIR "${MKLROOT_PATH}/compiler/lib/ia32")
        set(MKL_COMPILER_LIB_DIR ${MKL_COMPILER_LIB_DIR} "${MKLROOT_PATH}/lib/ia32")
        if ( WIN32 )
            set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_c)
        else ( WIN32 )
            set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel)
        endif ( WIN32 )
    endif ( USE_MKL_64BIT )

    if (WIN32)
        SET(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_rt)
        SET(MKL_LIBRARIES ${MKL_LIBRARIES} libiomp5md)
    else (WIN32)
        SET(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_thread)
        SET(MKL_LIBRARIES ${MKL_LIBRARIES} iomp5)
        SET(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_core)
    endif (WIN32) 
endif (MKL_FOUND)

IF (MKL_FOUND)
    IF (NOT MKL_FIND_QUIETLY)
        MESSAGE(STATUS "Found MKL libraries: ${MKL_LIBRARIES}")
        MESSAGE(STATUS "MKL_INCLUDE_DIR: ${MKL_INCLUDE_DIR}")
        MESSAGE(STATUS "MKL_LIB_DIR: ${MKL_LIB_DIR}")
        MESSAGE(STATUS "MKL_COMPILER_LIB_DIR: ${MKL_COMPILER_LIB_DIR}")
    ENDIF (NOT MKL_FIND_QUIETLY)

    # ------------------------------------------------------------------------
    #  Extract version information from <mkl.h>
    # ------------------------------------------------------------------------

    set(INTEL_MKL_VERSION_MAJOR 0)
    set(INTEL_MKL_VERSION_MINOR 0)
    set(INTEL_MKL_VERSION_UPDATE 0)

    if(EXISTS "${MKL_INCLUDE_DIR}/mkl.h")

        # Read and parse header file for version number
        file(STRINGS "${MKL_INCLUDE_DIR}/mkl.h" _mkl_HEADER_CONTENTS REGEX "#define __INTEL_MKL__ ")
        string(REGEX REPLACE ".*#define __INTEL_MKL__ ([0-9]+).*" "\\1" INTEL_MKL_VERSION_MAJOR "${_mkl_HEADER_CONTENTS}")
        unset(_mkl_HEADER_CONTENTS)
        
        file(STRINGS "${MKL_INCLUDE_DIR}/mkl.h" _mkl_HEADER_CONTENTS REGEX "#define __INTEL_MKL_MINOR__ ")
        string(REGEX REPLACE ".*#define __INTEL_MKL_MINOR__ ([0-9]+).*" "\\1" INTEL_MKL_VERSION_MINOR "${_mkl_HEADER_CONTENTS}")
        unset(_mkl_HEADER_CONTENTS)
        
        file(STRINGS "${MKL_INCLUDE_DIR}/mkl.h" _mkl_HEADER_CONTENTS REGEX "#define __INTEL_MKL_UPDATE__ ")
        string(REGEX REPLACE ".*#define __INTEL_MKL_UPDATE__ ([0-9]+).*" "\\1" INTEL_MKL_VERSION_UPDATE "${_mkl_HEADER_CONTENTS}")

        unset(_mkl_HEADER_CONTENTS)
    else (EXISTS "${MKL_INCLUDE_DIR}/mkl.h")
        if(EXISTS "${MKL_INCLUDE_DIR}/mkl_version.h")

            # Read and parse header file for version number
            file(STRINGS "${MKL_INCLUDE_DIR}/mkl_version.h" _mkl_HEADER_CONTENTS REGEX "#define __INTEL_MKL__ ")
            string(REGEX REPLACE ".*#define __INTEL_MKL__ ([0-9]+).*" "\\1" INTEL_MKL_VERSION_MAJOR "${_mkl_HEADER_CONTENTS}")
            unset(_mkl_HEADER_CONTENTS)
            
            file(STRINGS "${MKL_INCLUDE_DIR}/mkl_version.h" _mkl_HEADER_CONTENTS REGEX "#define __INTEL_MKL_MINOR__ ")
            string(REGEX REPLACE ".*#define __INTEL_MKL_MINOR__ ([0-9]+).*" "\\1" INTEL_MKL_VERSION_MINOR "${_mkl_HEADER_CONTENTS}")
            unset(_mkl_HEADER_CONTENTS)
            
            file(STRINGS "${MKL_INCLUDE_DIR}/mkl_version.h" _mkl_HEADER_CONTENTS REGEX "#define __INTEL_MKL_UPDATE__ ")
            string(REGEX REPLACE ".*#define __INTEL_MKL_UPDATE__ ([0-9]+).*" "\\1" INTEL_MKL_VERSION_UPDATE "${_mkl_HEADER_CONTENTS}")

            unset(_mkl_HEADER_CONTENTS)
        endif (EXISTS "${MKL_INCLUDE_DIR}/mkl_version.h")
    endif (EXISTS "${MKL_INCLUDE_DIR}/mkl.h")

    set(MKL_VERSION_STRING "${INTEL_MKL_VERSION_MAJOR}.${INTEL_MKL_VERSION_MINOR}.${INTEL_MKL_VERSION_UPDATE}")
    message("find MKL version : ${MKL_VERSION_STRING}")

    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
ELSE (MKL_FOUND)
    IF (MKL_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find MKL libraries")
    ENDIF (MKL_FIND_REQUIRED)
ENDIF (MKL_FOUND)

# MARK_AS_ADVANCED(MKL_LIBRARY)
