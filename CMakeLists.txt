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
if (OPENMP_FOUND)
  message("OpenMP multithreading enabled")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  OPTION(USE_OPENMP "Use OpenMP" On)
else (OPENMP_FOUND)
  message("OpenMP multithreading not supported")
endif (OPENMP_FOUND)
if(USE_OPENMP) 
  ADD_DEFINITIONS(-DUSE_OMP)
endif(USE_OPENMP) 

if (WIN32)
  ADD_DEFINITIONS(-DWIN32 -D_WIN32 -D_WINDOWS)
  ADD_DEFINITIONS(-DUNICODE -D_UNICODE)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
endif (WIN32)

set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/gadgetron)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

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

find_package(MKL)

if (CUDA_FOUND)
  MESSAGE("CUDA found, GPU components will be compiled.")
  SET( GADGETRON_CUDA_FOUND_BOOL 1 )
  include_directories( ${CUDA_INCLUDE_DIRS} )
  #set(CUDA_VERBOSE_BUILD ON)

  # Compile kernels for compute models 1.0 and 2.0 as default for Cuda 4.1
  # Support compute model 3.0 from Cuda 4.2 and up
  # Support compute model 3.5 from Cuda 5 and up

  set(CUDA_NVCC_FLAGS1 "-gencode arch=compute_10,code=sm_10")
  set(CUDA_NVCC_FLAGS2 "-gencode arch=compute_20,code=sm_20")
  set(CUDA_NVCC_FLAGS3 "-gencode arch=compute_30,code=sm_30") 
  set(CUDA_NVCC_FLAGS4 "-gencode arch=compute_35,code=sm_35")   

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
find_package(Numpy)
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

add_subdirectory(toolboxes)
add_subdirectory(apps)
if (ACE_FOUND AND ISMRMRD_FOUND)
  add_subdirectory(gadgets)
endif (ACE_FOUND AND ISMRMRD_FOUND)

add_subdirectory(cmake)
add_subdirectory(doc)

if (GTEST_FOUND AND CUDA_FOUND AND ARMADILLO_FOUND)
  add_subdirectory(test)
endif (GTEST_FOUND AND CUDA_FOUND AND ARMADILLO_FOUND)
