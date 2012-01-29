if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_HDF5MRI__)
endif (WIN32)

include_directories(${ACE_INCLUDE_DIR} 
				    ${Boost_INCLUDE_DIR}  
					${CMAKE_SOURCE_DIR}/toolboxes/ndarray
					${CMAKE_SOURCE_DIR}/toolboxes/hostutils
					${CMAKE_SOURCE_DIR}/toolboxes/hdf5utils
					${CMAKE_SOURCE_DIR}/gadgets/core)

add_library(hdf5mri SHARED mri_hdf5_io.cpp)

target_link_libraries(hdf5mri hdf5utils hondarray ${HDF5_LIBRARIES})


install(TARGETS hdf5mri DESTINATION lib)

install(FILES mri_hdf5_io.h hdf5mri_export.h DESTINATION include)
