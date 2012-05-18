if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON__)
ADD_DEFINITIONS(-DHDF5CPP_USEDLL)
ADD_DEFINITIONS(-D_HDF5USEDLL_)
endif (WIN32)

include_directories(${ACE_INCLUDE_DIR} 
				    ${Boost_INCLUDE_DIR}
					${HDF5_INCLUDE_DIRS}
					${HDF5_INCLUDE_DIRS}/cpp
					${CMAKE_SOURCE_DIR}/toolboxes/cpucore
					${CMAKE_SOURCE_DIR}/toolboxes/hostutils)

add_library(hdf5utils SHARED hdf5_core.cpp)
add_executable(hdf5_test hdf5_test.cpp)

target_link_libraries(hdf5utils cpucore ${HDF5_LIBRARIES})
target_link_libraries(hdf5_test hdf5utils cpucore ${HDF5_LIBRARIES})

if (WIN32)
target_link_libraries(hdf5utils ${HDF5_LIB_DIR}/libszip.lib ${HDF5_LIB_DIR}/libzlib.lib)
target_link_libraries(hdf5_test ${HDF5_LIB_DIR}/libszip.lib ${HDF5_LIB_DIR}/libzlib.lib)
endif (WIN32)

install(TARGETS hdf5utils DESTINATION lib)
install(TARGETS hdf5_test DESTINATION bin)

install(FILES hdf5_core.h hoNDArray_hdf5_io.h hdf5utils_export.h DESTINATION include)
