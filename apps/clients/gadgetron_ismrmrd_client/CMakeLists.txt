find_package(ZFP)


if (ZFP_FOUND)
   add_definitions(-DGADGETRON_COMPRESSION_ZFP)
endif()

if (ZFP_FOUND)
   include_directories(${ZFP_INCLUDE_DIR})
endif()


add_executable(gadgetron_ismrmrd_client gadgetron_ismrmrd_client.cpp)
include_directories(${HDF5_INCLUDE_DIR})
target_link_libraries(gadgetron_ismrmrd_client ISMRMRD::ISMRMRD boost gadgetron_mricore )

if (ZFP_FOUND)
   target_link_libraries(gadgetron_ismrmrd_client ${ZFP_LIBRARIES})
endif ()

install(TARGETS gadgetron_ismrmrd_client DESTINATION bin COMPONENT main)
