cmake_minimum_required(VERSION 2.6)

set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/gadgetron)

#set(CMAKE_OSX_ARCHITECTURES "i386")

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

project(GADGETRON)

find_package(Boost REQUIRED)
find_package(FFTW3 COMPONENTS single double REQUIRED)
find_package(ACE REQUIRED)
find_package(TinyXML REQUIRED)

find_package(CUDA 4.0)

if (CUDA_FOUND)
   if (CUDA_SDK_ROOT_DIR MATCHES "CUDA_SDK_ROOT_DIR-NOTFOUND")
      MESSAGE(FATAL_ERROR "Could not locate CUDA SDK. Try setting environment variable NVSDKCOMPUTE_ROOT")
    endif (CUDA_SDK_ROOT_DIR MATCHES "CUDA_SDK_ROOT_DIR-NOTFOUND")
else (CUDA_FOUND)
     MESSAGE("CUDA not found. CUDA components will not be compiled")
endif (CUDA_FOUND)

find_package(Qt4 4.7)

find_package(PythonLibs)

set(CMAKE_CXX_FLAGS "-g -Wall")

add_subdirectory(apps)
add_subdirectory(toolboxes)
add_subdirectory(gadgets)
