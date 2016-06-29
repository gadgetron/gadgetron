if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_INTERVENTIONAL_MRI__)
endif (WIN32)

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
    ${HDF5_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}/cpp
    ${ARMADILLO_INCLUDE_DIRS}
    ${MKL_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
)


add_library(gadgetron_interventional_mri SHARED 
    gadgetron_interventional_mri_export.h 
    DeviceChannelSplitterGadget.h
    DeviceChannelSplitterGadget.cpp
)

set_target_properties(gadgetron_interventional_mri PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})                                                                                                                                                                                                      

target_link_libraries(gadgetron_interventional_mri
    gadgetron_gadgetbase
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_gadgettools
    ${ISMRMRD_LIBRARIES} 
    ${ACE_LIBRARIES}
)

install(FILES 
    gadgetron_interventional_mri_export.h
    DeviceChannelSplitterGadget.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES grappa_device.xml grappa_device_cpu.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install(TARGETS gadgetron_interventional_mri DESTINATION lib COMPONENT main)
