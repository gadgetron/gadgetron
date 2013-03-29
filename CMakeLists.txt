cmake_minimum_required(VERSION 2.8)
project(GADGETRON)

if (WIN32)
  ADD_DEFINITIONS(-DWIN32 -D_WIN32 -D_WINDOWS)
  ADD_DEFINITIONS(-DUNICODE -D_UNICODE)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
endif (WIN32)

set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/gadgetron)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED)
find_package(ACE REQUIRED)
find_package(FFTW3 COMPONENTS single double REQUIRED)

find_package(CUDA 4.1)
find_package(GTest)
find_package(Armadillo)

if (CUDA_FOUND)
  MESSAGE("CUDA Found, GPU components will be compiled")
  SET( GADGETRON_CUDA_FOUND_BOOL 1 )
  include_directories( ${CUDA_INCLUDE_DIRS} )
else (CUDA_FOUND)
  MESSAGE("CUDA not found. CUDA components will not be compiled")
  SET( GADGETRON_CUDA_FOUND_BOOL 0 )
endif (CUDA_FOUND)

# Compile kernels for shader models 1.0 and 2.0 as default
IF (${CMAKE_VERSION} VERSION_GREATER "2.8.7")
  set(CUDA_NVCC_FLAGS "-gencode=arch=compute_10,code=\"sm_10,compute_10\" -gencode=arch=compute_20,code=\"sm_20,compute_20\"")
ENDIF (${CMAKE_VERSION} VERSION_GREATER "2.8.7")

find_package(Qt4 4.6)

find_package(PythonLibs)
find_package(Numpy)
find_package(GLEW)
find_package(Ismrmrd)

add_subdirectory(dependencies)
add_subdirectory(toolboxes)
add_subdirectory(apps)
IF (ACE_FOUND AND ISMRMRD_FOUND)
  add_subdirectory(gadgets)
ENDIF (ACE_FOUND AND ISMRMRD_FOUND)

add_subdirectory(cmake)
add_subdirectory(doc)

if (GTEST_FOUND)
  ENABLE_TESTING()
  add_subdirectory(test)
endif (GTEST_FOUND)
