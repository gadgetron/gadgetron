if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_INTERVENTIONAL_MRI__)
endif ()



add_library(gadgetron_interventional_mri SHARED 
    gadgetron_interventional_mri_export.h 
    DeviceChannelSplitterGadget.h
    DeviceChannelSplitterGadget.cpp
)

set_target_properties(gadgetron_interventional_mri PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})                                                                                                                                                                                                      

target_link_libraries(gadgetron_interventional_mri
    gadgetron_core
    gadgetron_mricore
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
)

install(FILES
    gadgetron_interventional_mri_export.h
    DeviceChannelSplitterGadget.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES grappa_device.xml grappa_device_cpu.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install(TARGETS gadgetron_interventional_mri
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)