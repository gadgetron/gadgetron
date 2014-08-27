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
  ${ACE_INCLUDE_DIR} 
  ${Boost_INCLUDE_DIR} 
  ${HDF5_CXX_INCLUDE_DIR} 
  ${HDF5_C_INCLUDE_DIR} 
  ${ISMRMRD_INCLUDE_DIR}
  )

add_executable(ismrmrd_client ismrmrd_client.cpp)

target_link_libraries(ismrmrd_client ${HDF5_LIBRARIES} ${ISMRMRD_LIBRARIES} ${Boost_LIBRARIES})

install(TARGETS ismrmrd_client DESTINATION bin)
