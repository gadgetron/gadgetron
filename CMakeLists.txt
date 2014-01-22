cmake_minimum_required(VERSION 2.8)
project(GADGETRON)

# build options for 64 bits system
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  message(" 64bit system is found  ... ")
  set( HAS_64_BIT On CACHE BOOL "64bit build")
else( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  message(" 32bit system is found  ... ")
  set( HAS_64_BIT Off CACHE BOOL "64bit build"e)
endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )

# build options for OpenMP support
find_package(OpenMP)
OPTION(USE_OPENMP "Use OpenMP" On)
if (OPENMP_FOUND)
    if(USE_OPENMP) 
        message("OpenMP multithreading enabled")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        ADD_DEFINITIONS(-DUSE_OMP)
    else (USE_OPENMP)
        message("OpenMP multithreading is supported, but disabled")
    endif(USE_OPENMP) 
else (OPENMP_FOUND)
  message("OpenMP multithreading not supported")
endif (OPENMP_FOUND)

if (WIN32)
    ADD_DEFINITIONS(-DWIN32 -D_WIN32 -D_WINDOWS)
    ADD_DEFINITIONS(-DUNICODE -D_UNICODE)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    if ( HAS_64_BIT )
        ADD_DEFINITIONS(-DWIN64 -D_WIN64)
    endif ( HAS_64_BIT )
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
else (WIN32)
    if (UNIX)
        if (APPLE)
        else (APPLE)
            SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
            SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libgfortran")
        endif (APPLE)
  endif (UNIX)
endif (WIN32)

set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/gadgetron)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/prod/bin/vc10 CACHE STRING "Where to put the executables")
SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/prod/bin/vc10 CACHE STRING "Where to put the libraries")

# whether to suppress compilation warnings
OPTION(BUILD_SUPPRESS_WARNINGS "Build package while suppressing warnings" Off)
if (BUILD_SUPPRESS_WARNINGS)
    if (WIN32)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    elseif (WIN32)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w")
    endif (WIN32)
endif (BUILD_SUPPRESS_WARNINGS)

# whether to compile toolboxes as static library
OPTION(BUILD_TOOLBOX_STATIC "Build static library for toolboxes" Off)

if ( BUILD_TOOLBOX_STATIC )
    message("Build static toolbox libray ... ")
    ADD_DEFINITIONS(-DBUILD_TOOLBOX_STATIC)
    set(LIBTYPE STATIC)
    if ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
    endif ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
else ( BUILD_TOOLBOX_STATIC )
    message("Build dynamic toolbox libray ... ")
    set(LIBTYPE SHARED)
endif ( BUILD_TOOLBOX_STATIC )

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_NO_BOOST_CMAKE ON)

# We actually only use system and thread explicitly, but they require linking in date_time and chrono...
if (WIN32)
  find_package(Boost COMPONENTS system thread date_time chrono REQUIRED)
else(WIN32)
  find_package(Boost COMPONENTS system thread REQUIRED)
endif(WIN32)

find_package(FFTW3 COMPONENTS single double REQUIRED)

find_package(ACE)
if(ACE_FOUND)
  MESSAGE("ACE found, the streaming framework will be compiled.")
else(ACE_FOUND)
  MESSAGE("ACE not found. Only toolboxes and standalone applications are compiled. The streaming framework will not be compiled.")
endif(ACE_FOUND)

find_package(CUDA 4.1)
if ( CUDA_FOUND )
    ADD_DEFINITIONS(-DUSE_CUDA)
endif ( CUDA_FOUND )

find_package(GTest)
#Add support for the default ubuntu package of gtest (which is not compiled
if (NOT GTEST_FOUND)
  find_path(GTEST_SRC_DIR src/gtest.cc HINTS /usr/src/gtest)
  find_path(GTEST_INCLUDE_DIRS gtest.h HINTS /usr/include/gtest)
  if (GTEST_SRC_DIR AND GTEST_INCLUDE_DIRS)
    MESSAGE("GTest src package found. Compiling as part of Gadgetron.")
    add_subdirectory(${GTEST_SRC_DIR} ${CMAKE_BINARY_DIR}/gtest )
    include_directories(${GTEST_INCLUDE_DIRS})
    set(GTEST_FOUND 1)
    set(GTEST_LIBRARIES gtest gtest_main)
  endif (GTEST_SRC_DIR AND GTEST_INCLUDE_DIRS)
endif (NOT GTEST_FOUND)

find_package(Armadillo)
# check whether ILP64 MKL should is used
if(ARMADILLO_FOUND)
    set(ARMADILLO_BLAS_LONG_LONG FALSE)
    if(EXISTS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp")
        # Read and parse armadillo config.hpp to find out whether BLAS uses long long
        file(STRINGS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp" _armadillo_blas_long_long REGEX "// #define ARMA_BLAS_LONG_LON")
        if ( NOT _armadillo_blas_long_long )
            set(ARMADILLO_BLAS_LONG_LONG TRUE)
            MESSAGE("Armadillo is found to use long long for BLAS calls ... ")
        else ( NOT _armadillo_blas_long_long )
            MESSAGE("Armadillo is found NOT to use long long for BLAS calls ... ")
            MESSAGE("Note the ARMADILLO_BLAS_LONG_LONG must be defined in the ${ARMADILLO_INCLUDE_DIR}/armadillo_bits/config.hpp to link against MKL ILP64 interface ... ")
        endif ( NOT _armadillo_blas_long_long )
        unset(_armadillo_blas_long_long)
    endif()
endif ()

find_package(MKL)

if (CUDA_FOUND)
  MESSAGE("CUDA found, GPU components will be compiled.")
  SET( GADGETRON_CUDA_FOUND_BOOL 1 )
  include_directories( ${CUDA_INCLUDE_DIRS} )
  #set(CUDA_VERBOSE_BUILD ON)

  # Compile kernels for compute models 1.0 and 2.0 as default for Cuda 4.1
  # Support compute model 3.0 from Cuda 4.2 and up
  # Support compute model 3.5 from Cuda 5 and up

  #set(CUDA_NVCC_FLAGS1 "-gencode arch=compute_10,code=sm_10")
  #set(CUDA_NVCC_FLAGS2 "-gencode arch=compute_20,code=sm_20")
  set(CUDA_NVCC_FLAGS3 "-gencode arch=compute_30,code=sm_30") 
  #set(CUDA_NVCC_FLAGS4 "-gencode arch=compute_35,code=sm_35")   

  if(${CUDA_VERSION} VERSION_GREATER "4.99")
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS1} ${CUDA_NVCC_FLAGS2} ${CUDA_NVCC_FLAGS3} ${CUDA_NVCC_FLAGS4})
  else(${CUDA_VERSION} VERSION_GREATER "4.99")    

    if(${CUDA_VERSION} VERSION_GREATER "4.1")
      set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS1} ${CUDA_NVCC_FLAGS2} ${CUDA_NVCC_FLAGS3})
    else(${CUDA_VERSION} VERSION_GREATER "4.1")      

      set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS1} ${CUDA_NVCC_FLAGS2})          

    endif(${CUDA_VERSION} VERSION_GREATER "4.1")
  endif(${CUDA_VERSION} VERSION_GREATER "4.99")

else (CUDA_FOUND)
  MESSAGE("CUDA not found. CUDA components will not be compiled.")
  SET( GADGETRON_CUDA_FOUND_BOOL 0 )
endif (CUDA_FOUND)

find_package(Qt4 4.6)
find_package(PythonLibs)
find_package(NumPy)
find_package(GLEW)
find_package(OpenGL)
find_package(GLUT)
find_package(HDF5 1.8 COMPONENTS C CXX)

find_package(Ismrmrd)
if(ISMRMRD_FOUND)
  message("ISMRMRD found")
  find_package(XSD REQUIRED)
  find_package(XercesC REQUIRED)
else(ISMRMRD_FOUND)
  message("ISMRMRD not found. Only compiling toolboxes and standalone applications.")
endif(ISMRMRD_FOUND)

find_package(GMatlab)

include_directories( ${CMAKE_SOURCE_DIR} )

add_subdirectory(toolboxes)
add_subdirectory(apps)
if (ACE_FOUND AND ISMRMRD_FOUND)
  add_subdirectory(gadgets)
endif (ACE_FOUND AND ISMRMRD_FOUND)

add_subdirectory(cmake)
add_subdirectory(doc)

if (GTEST_FOUND AND ARMADILLO_FOUND)
  add_subdirectory(test)
endif (GTEST_FOUND AND ARMADILLO_FOUND)
