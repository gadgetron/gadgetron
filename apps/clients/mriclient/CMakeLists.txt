find_package(Ismrmrd REQUIRED)

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
  ${CMAKE_SOURCE_DIR}/toolboxes/hdf5utils
  ${CMAKE_SOURCE_DIR}/toolboxes/hdf5mri
  ${ISMRMRD_INCLUDE_DIR}
  )

add_executable(mriclient main.cpp)

target_link_libraries(mriclient tinyxml gadgettools hdf5utils hdf5mri optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${ISMRMRD_LIBRARIES})

if (WIN32)
    target_link_libraries(mriclient optimized ${HDF5_LIB_DIR}/hdf5dll.lib)
    target_link_libraries(mriclient optimized ${HDF5_LIB_DIR}/hdf5_cppdll.lib)
    target_link_libraries(mriclient optimized ${HDF5_LIB_DIR}/hdf5_hldll.lib)

    target_link_libraries(mriclient debug ${HDF5_LIB_DIR}/hdf5ddll.lib)
    target_link_libraries(mriclient debug ${HDF5_LIB_DIR}/hdf5_cppddll.lib)
    target_link_libraries(mriclient debug ${HDF5_LIB_DIR}/hdf5_hlddll.lib)


    target_link_libraries(mriclient optimized ${HDF5_LIB_DIR}/zlib.lib)
    target_link_libraries(mriclient optimized ${HDF5_LIB_DIR}/szip.lib)
    target_link_libraries(mriclient debug ${HDF5_LIB_DIR}/zlib_D.lib)
    target_link_libraries(mriclient debug ${HDF5_LIB_DIR}/szip_D.lib)

else (WIN32)

    target_link_libraries(mriclient ${HDF5_LIBRARIES})

endif (WIN32)

install(TARGETS mriclient DESTINATION bin)
install (FILES ImageWriter.h HDF5ImageWriter.h DESTINATION include)
install(FILES ${ISMRMRD_LIBRARIES} DESTINATION lib)
