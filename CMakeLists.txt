cmake_minimum_required(VERSION 2.6)

PROJECT(gt-tomography)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
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

