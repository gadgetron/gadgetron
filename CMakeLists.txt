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
   MESSAGE("CUDA Found, GPU components will be compiled")
else (CUDA_FOUND)
     MESSAGE("CUDA not found. CUDA components will not be compiled")
endif (CUDA_FOUND)

find_package(Qt4 4.7)

find_package(PythonLibs)


find_package(Doxygen)
	if(DOXYGEN_FOUND)
		configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
		add_custom_target(doc ${DOXYGEN_EXECUTABLE}
	${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM)
else(DOXYGEN_FOUND)
	MESSAGE("Doxygen not found. Will not be able to build documentation")
endif(DOXYGEN_FOUND)

set(CMAKE_CXX_FLAGS "-g -Wall")

add_subdirectory(apps)
add_subdirectory(toolboxes)
add_subdirectory(gadgets)
