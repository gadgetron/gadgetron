find_package(Ismrmrd REQUIRED)
find_package(HDF5 1.8 COMPONENTS C CXX REQUIRED)
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS thread system REQUIRED)

include_directories(    
  ${ACE_INCLUDE_DIR} 
  ${Boost_INCLUDE_DIR}  
  ${HDF5_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}/cpp
  ${CMAKE_SOURCE_DIR}/gadgets/core 
  ${CMAKE_SOURCE_DIR}/apps/gadgetron
  ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
  ${CMAKE_SOURCE_DIR}/toolboxes/hostutils
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  ${ISMRMRD_INCLUDE_DIR}
  ${XSD_INCLUDE_DIR}
  )

LINK_DIRECTORIES(${BOOST_ROOT}/lib)

add_executable(mriclient main.cpp)

IF(WIN32)
	SET(HDF5_LIB_DIR ${HDF5_INCLUDE_DIR}/../lib)
	target_link_libraries(mriclient optimized ${HDF5_LIB_DIR}/hdf5dll.lib)
	target_link_libraries(mriclient optimized ${HDF5_LIB_DIR}/hdf5_cppdll.lib)
	target_link_libraries(mriclient optimized ${HDF5_LIB_DIR}/hdf5_hldll.lib)

	target_link_libraries(mriclient debug ${HDF5_LIB_DIR}/hdf5ddll.lib)
	target_link_libraries(mriclient debug ${HDF5_LIB_DIR}/hdf5_cppddll.lib)
	target_link_libraries(mriclient debug ${HDF5_LIB_DIR}/hdf5_hlddll.lib)

	#target_link_libraries(mriclient optimized ${HDF5_LIB_DIR}/zlib.lib)
	#target_link_libraries(mriclient optimized ${HDF5_LIB_DIR}/szip.lib)
	#target_link_libraries(mriclient debug ${HDF5_LIB_DIR}/zlib.lib)
	#target_link_libraries(mriclient debug ${HDF5_LIB_DIR}/szip.lib)

	target_link_libraries(mriclient tinyxml gadgettools optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${ISMRMRD_LIBRARIES}  ${ISMRMRD_XSD_LIBRARIES} ${Boost_LIBRARIES})
ELSE (WIN32)
	target_link_libraries(mriclient tinyxml gadgettools  ${HDF5_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${ISMRMRD_LIBRARIES}  ${ISMRMRD_XSD_LIBRARIES} ${Boost_LIBRARIES})
ENDIF(WIN32)


install(TARGETS mriclient DESTINATION bin)
install (FILES ImageWriter.h HDF5ImageWriter.h BlobFileWriter.h DESTINATION include)
install(FILES ${ISMRMRD_LIBRARIES} DESTINATION lib)
