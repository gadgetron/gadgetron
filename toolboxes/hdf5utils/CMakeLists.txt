if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON__)
ADD_DEFINITIONS(-DHDF5CPP_USEDLL)
ADD_DEFINITIONS(-D_HDF5USEDLL_)
ADD_DEFINITIONS(-D_HDF5USEHLDLL_)
ADD_DEFINITIONS(-DHDF5USE_HLCPPDLL)
endif (WIN32)

message(${HDF5_LIBRARIES})

include_directories(${ACE_INCLUDE_DIR} 
				    ${Boost_INCLUDE_DIR}
					${HDF5_INCLUDE_DIRS}
					${HDF5_INCLUDE_DIRS}/cpp
					${CMAKE_SOURCE_DIR}/toolboxes/cpucore
					${CMAKE_SOURCE_DIR}/toolboxes/hostutils)

add_library(hdf5utils SHARED hdf5_core.cpp)
add_executable(hdf5_test hdf5_test.cpp)

# target_link_libraries(hdf5utils cpucore ${HDF5_LIBRARIES})
# target_link_libraries(hdf5_test hdf5utils cpucore ${HDF5_LIBRARIES})

#if (WIN32)
#target_link_libraries(hdf5utils ${HDF5_LIB_DIR}/zlib.lib ${HDF5_LIB_DIR}/szip.lib)
#target_link_libraries(hdf5_test ${HDF5_LIB_DIR}/zlib.lib ${HDF5_LIB_DIR}/szip.lib)
#endif (WIN32)

if (WIN32)

target_link_libraries(hdf5utils cpucore gadgetroncore)
target_link_libraries(hdf5_test hdf5utils cpucore gadgetroncore)

target_link_libraries(hdf5utils optimized ${HDF5_LIB_DIR}/hdf5dll.lib)
target_link_libraries(hdf5utils optimized ${HDF5_LIB_DIR}/hdf5_cppdll.lib)
target_link_libraries(hdf5utils optimized ${HDF5_LIB_DIR}/hdf5_hldll.lib)

target_link_libraries(hdf5utils debug ${HDF5_LIB_DIR}/hdf5ddll.lib)
target_link_libraries(hdf5utils debug ${HDF5_LIB_DIR}/hdf5_cppddll.lib)
target_link_libraries(hdf5utils debug ${HDF5_LIB_DIR}/hdf5_hlddll.lib)

target_link_libraries(hdf5utils optimized ${HDF5_LIB_DIR}/zlib.lib)
target_link_libraries(hdf5utils optimized ${HDF5_LIB_DIR}/szip.lib)
target_link_libraries(hdf5utils debug ${HDF5_LIB_DIR}/zlib_D.lib)
target_link_libraries(hdf5utils debug ${HDF5_LIB_DIR}/szip_D.lib)

target_link_libraries(hdf5_test optimized ${HDF5_LIB_DIR}/zlib.lib)
target_link_libraries(hdf5_test optimized ${HDF5_LIB_DIR}/szip.lib)
target_link_libraries(hdf5_test debug ${HDF5_LIB_DIR}/zlib_D.lib)
target_link_libraries(hdf5_test debug ${HDF5_LIB_DIR}/szip_D.lib)

else (WIN32)

target_link_libraries(hdf5utils cpucore ${HDF5_LIBRARIES})
target_link_libraries(hdf5_test hdf5utils cpucore ${HDF5_LIBRARIES})

endif (WIN32)

install(TARGETS hdf5utils DESTINATION lib)
install(TARGETS hdf5_test DESTINATION bin)

install(FILES hdf5_core.h hoNDArray_hdf5_io.h hdf5utils_export.h DESTINATION include)
