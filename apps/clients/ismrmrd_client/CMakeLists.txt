find_package(Ismrmrd REQUIRED)
find_package(HDF5 1.8 COMPONENTS C CXX REQUIRED)

set(Boost_NO_BOOST_CMAKE ON)

if(WIN32)
  find_package(Boost COMPONENTS program_options thread system date_time chrono REQUIRED)
else(WIN32)
  find_package(Boost COMPONENTS program_options thread system REQUIRED)
endif(WIN32)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
  ${ISMRMRD_SCHEMA_DIR}
  ${ISMRMRD_XSD_INCLUDE_DIR}
  ${ACE_INCLUDE_DIR} 
  ${Boost_INCLUDE_DIR} 
  ${HDF5_CXX_INCLUDE_DIR} 
  ${HDF5_C_INCLUDE_DIR} 
  ${ISMRMRD_INCLUDE_DIR}
  )

add_executable(ismrmrd_client ismrmrd_client.cpp)

IF(WIN32)
    target_link_libraries(ismrmrd_client optimized ${HDF5_hdf5_LIBRARY_RELEASE} optimized ${HDF5_hdf5_cpp_LIBRARY_RELEASE})
    target_link_libraries(ismrmrd_client debug ${HDF5_hdf5_LIBRARY_DEBUG} debug ${HDF5_hdf5_cpp_LIBRARY_DEBUG})
    target_link_libraries(ismrmrd_client ${ISMRMRD_LIBRARIES} ${Boost_LIBRARIES})
ELSE (WIN32)
    target_link_libraries(ismrmrd_client ${HDF5_LIBRARIES} ${ISMRMRD_LIBRARIES} ${Boost_LIBRARIES})
ENDIF(WIN32)


install(TARGETS ismrmrd_client DESTINATION bin)
