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
find_package(Boost REQUIRED)
find_package(CUDA 5.0)
if (CUDA_FOUND)
   MESSAGE("CUDA Found, GPU components will be compiled")
   ADD_DEFINITIONS(-DGTCUDA_FOUND)
   include_directories( ${CUDA_INCLUDE_DIRS})
   set(CUDA_NVCC_FLAGS "-arch=sm_20")
else (CUDA_FOUND)
     MESSAGE("CUDA not found. CUDA components will not be compiled")
endif (CUDA_FOUND)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
add_subdirectory(solvers)
add_subdirectory(operators)
add_subdirectory(proton)

