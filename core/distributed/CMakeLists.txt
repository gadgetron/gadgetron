
add_library(gadgetron_core_distributed SHARED
        Distributor.h
        ChannelCreator.h AcquisitionDistributor.cpp AcquisitionDistributor.h Distributor.cpp BufferDistributor.cpp BufferDistributor.h)


target_link_libraries(gadgetron_core_distributed
        gadgetron_core
        ismrmrd
        gadgetron_toolbox_mri_core)

set_target_properties(gadgetron_core_distributed PROPERTIES
        VERSION ${GADGETRON_VERSION_STRING}
        SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS gadgetron_core_distributed
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)