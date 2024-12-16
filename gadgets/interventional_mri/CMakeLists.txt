add_library(pingvin_interventional_mri SHARED
    DeviceChannelSplitterGadget.h
    DeviceChannelSplitterGadget.cpp
)

set_target_properties(pingvin_interventional_mri PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_interventional_mri
    pingvin_core
    pingvin_mricore
    pingvin_toolbox_log
    pingvin_toolbox_cpucore
)

install(FILES
    DeviceChannelSplitterGadget.h
    DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES grappa_device.xml grappa_device_cpu.xml DESTINATION ${PINGVIN_INSTALL_CONFIG_PATH} COMPONENT main)

install(TARGETS pingvin_interventional_mri
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)