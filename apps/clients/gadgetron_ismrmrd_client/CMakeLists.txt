find_package(ZFP)

if (ZFP_FOUND)
   add_definitions(-DGADGETRON_COMPRESSION)
endif(ZFP_FOUND)

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
  ${Boost_INCLUDE_DIR} 
  ${ISMRMRD_INCLUDE_DIR}
  )

if (ZFP_FOUND)
   include_directories(${ZFP_INCLUDE_DIR})
endif(ZFP_FOUND)

add_executable(gadgetron_ismrmrd_client gadgetron_ismrmrd_client.cpp)

target_link_libraries(gadgetron_ismrmrd_client ${ISMRMRD_LIBRARIES} ${Boost_LIBRARIES})

if (ZFP_FOUND)
   target_link_libraries(gadgetron_ismrmrd_client ${ZFP_LIBRARIES})
endif (ZFP_FOUND)

install(TARGETS gadgetron_ismrmrd_client DESTINATION bin COMPONENT main)
