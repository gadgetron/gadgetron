cmake_minimum_required(VERSION 2.6)

PROJECT(gt-tomography)


find_path(GADGETRON_CMAKE_MODULES FindGadgetron.cmake HINTS
$ENV{GADGETRON_HOME}/cmake
/usr/local/gadgetron/cmake)

if (NOT GADGETRON_CMAKE_MODULES)
  MESSAGE(FATAL_ERROR "GADGETRON_CMAKE_MODULES cannot be found. 
   Try to set GADGETRON_HOME environment variable.")
endif(NOT GADGETRON_CMAKE_MODULES)

set(CMAKE_MODULE_PATH ${GADGETRON_CMAKE_MODULES})

find_package(Gadgetron REQUIRED)

include_directories(  "${GADGETRON_HOME}/include")

LINK_DIRECTORIES(${GADGETRON_LIB_DIR})   
find_package(Boost COMPONENTS program_options system thread REQUIRED)
find_package(CUDA 4.0 REQUIRED)
find_package(CULA REQUIRED)
if (CUDA_FOUND)
   MESSAGE("CUDA Found, GPU components will be compiled")
   ADD_DEFINITIONS(-DGTCUDA_FOUND)
   include_directories( ${CUDA_INCLUDE_DIRS})
   set(CUDA_NVCC_FLAGS "-arch=sm_20")
else (CUDA_FOUND)
     MESSAGE("CUDA not found. CUDA components will not be compiled")
endif (CUDA_FOUND)


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

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
add_subdirectory(solvers)
add_subdirectory(operators)
add_subdirectory(proton)
add_subdirectory(Xray)
add_subdirectory(test)
add_subdirectory(mri)

