if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_HDF5MRI__)
ADD_DEFINITIONS(-DHDF5CPP_USEDLL)
ADD_DEFINITIONS(-D_HDF5USEDLL_)
ADD_DEFINITIONS(-D_HDF5USEHLDLL_)
ADD_DEFINITIONS(-DHDF5USE_HLCPPDLL)
endif (WIN32)

include_directories(${ACE_INCLUDE_DIR} 
				    ${Boost_INCLUDE_DIR}
					${HDF5_INCLUDE_DIRS}
					${HDF5_INCLUDE_DIRS}/cpp
					${CMAKE_SOURCE_DIR}/toolboxes/cpucore
					${CMAKE_SOURCE_DIR}/toolboxes/hostutils
					${CMAKE_SOURCE_DIR}/toolboxes/hdf5utils
					${CMAKE_SOURCE_DIR}/gadgets/core)

add_library(hdf5mri SHARED mri_hdf5_io.cpp)

target_link_libraries(hdf5mri hdf5utils cpucore ${HDF5_LIBRARIES})

if (WIN32)
target_link_libraries(hdf5mri ${HDF5_LIB_DIR}/libszip.lib ${HDF5_LIB_DIR}/libzlib.lib)
endif (WIN32)

install(TARGETS hdf5mri DESTINATION lib)

install(FILES mri_hdf5_io.h hdf5mri_export.h DESTINATION include)
