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
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(Boost REQUIRED)
find_package(ACE REQUIRED)
find_package(FFTW3 COMPONENTS single double REQUIRED)

if (${CMAKE_VERSION} VERSION_LESS "2.8.4")
 MESSAGE("You have cmake version less than 2.8.4, hence FindBLAS requires a Fortran compiler")
 MESSAGE("Upgrade to a newer version of cmake to avoid this Fortran compiler check")
 enable_language(Fortran)
endif(${CMAKE_VERSION} VERSION_LESS "2.8.4")

find_package(BLAS)
find_package(LAPACK) 

find_package(CUDA 4.1)

if (CUDA_FOUND)
   MESSAGE("CUDA Found, GPU components will be compiled")
else (CUDA_FOUND)
     MESSAGE("CUDA not found. CUDA components will not be compiled")
endif (CUDA_FOUND)

find_package(Qt4 4.6)

find_package(PythonLibs)
find_package(Numpy)
find_package(Boost COMPONENTS python)

find_package(GLEW)

add_subdirectory(doc)
add_subdirectory(dependencies)
add_subdirectory(apps)
add_subdirectory(toolboxes)
IF (ACE_FOUND)
       add_subdirectory(gadgets)
ENDIF (ACE_FOUND)
add_subdirectory(cmake)
