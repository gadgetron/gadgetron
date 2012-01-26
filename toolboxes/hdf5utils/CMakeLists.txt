if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_HDF5UTILS__)
endif (WIN32)

include_directories(${ACE_INCLUDE_DIR} 
				    ${Boost_INCLUDE_DIR}  
					${CMAKE_SOURCE_DIR}/toolboxes/ndarray
					${CMAKE_SOURCE_DIR}/toolboxes/hostutils)

add_library(hdf5utils SHARED hoNDArray_hdf5_io.cpp)
add_executable(hdf5_test hdf5_test.cpp)

target_link_libraries(hdf5utils hondarray ${HDF5_LIBRARIES})
target_link_libraries(hdf5_test hdf5utils hondarray ${HDF5_LIBRARIES})


install(TARGETS hdf5utils DESTINATION lib)
install(TARGETS hdf5_test DESTINATION bin)

install(FILES hoNDArray_hdf5_io.h DESTINATION include)
