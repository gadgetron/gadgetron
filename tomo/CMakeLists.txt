cmake_minimum_required(VERSION 2.6)

PROJECT(gt-tomography)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(Gadgetron REQUIRED)



#find_path(GADGETRON_CMAKE_MODULES FindGadgetron.cmake HINTS
#$ENV{GADGETRON_HOME}/cmake
#/usr/local/gadgetron/cmake)

#if (NOT GADGETRON_CMAKE_MODULES)
#  MESSAGE(FATAL_ERROR "GADGETRON_CMAKE_MODULES cannot be found. 
#   Try to set GADGETRON_HOME environment variable.")
#endif(NOT GADGETRON_CMAKE_MODULES)


include_directories(  "${GADGETRON_HOME}/include")

LINK_DIRECTORIES(${GADGETRON_LIB_DIR})   
find_package(Boost COMPONENTS program_options system thread REQUIRED)
find_package(CUDA 4.0 REQUIRED)
find_package(CULA REQUIRED)
if (CUDA_FOUND)
   MESSAGE("CUDA Found, GPU components will be compiled")
   ADD_DEFINITIONS(-DGTCUDA_FOUND)
   include_directories( ${CUDA_INCLUDE_DIRS})
  # Compile kernels for shader models 2.0 for Cuda 4
  # Add shader models 3.0 and 3.5 from Cuda 5 and up
  set(CUDA_NVCC_FLAGS2 "-gencode arch=compute_20,code=sm_20")
  set(CUDA_NVCC_FLAGS3 "-gencode arch=compute_30,code=sm_30") 
  set(CUDA_NVCC_FLAGS4 "-gencode arch=compute_35,code=sm_35")   

  if(${CUDA_VERSION} VERSION_GREATER "4.99")
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS1} ${CUDA_NVCC_FLAGS2} ${CUDA_NVCC_FLAGS3} ${CUDA_NVCC_FLAGS4})    
  else(${CUDA_VERSION} VERSION_GREATER "4.99")
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS1} ${CUDA_NVCC_FLAGS2})    
  endif(${CUDA_VERSION} VERSION_GREATER "4.99")
else (CUDA_FOUND)
     MESSAGE("CUDA not found. CUDA components will not be compiled")
endif (CUDA_FOUND)

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
add_subdirectory(xray)
add_subdirectory(test)
add_subdirectory(mri)

