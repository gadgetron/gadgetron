find_package(ZFP)
find_package(HDF5 1.8 REQUIRED COMPONENTS C)

include_directories(${HDF5_C_INCLUDE_DIR})

if (ZFP_FOUND)
   add_definitions(-DGADGETRON_COMPRESSION_ZFP)
endif()

set(Boost_NO_BOOST_CMAKE ON)

if(WIN32)
  find_package(Boost COMPONENTS program_options thread system date_time chrono REQUIRED)
else()
  find_package(Boost COMPONENTS program_options thread system REQUIRED)
endif()

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif()

include_directories(
  ${Boost_INCLUDE_DIR} 
  ${ISMRMRD_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  )

if (ZFP_FOUND)
   include_directories(${ZFP_INCLUDE_DIR})
endif()


add_executable(gadgetron_ismrmrd_client gadgetron_ismrmrd_client.cpp)

target_link_libraries(gadgetron_ismrmrd_client ${ISMRMRD_LIBRARIES} ${Boost_LIBRARIES})

if (ZFP_FOUND)
   target_link_libraries(gadgetron_ismrmrd_client ${ZFP_LIBRARIES})
endif ()

install(TARGETS gadgetron_ismrmrd_client DESTINATION bin COMPONENT main)
