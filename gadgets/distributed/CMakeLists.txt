IF (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_DISTRIBUTED_GADGETS__)
ENDIF (WIN32)

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/cloudbus
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
)

add_library(gadgetron_distributed SHARED 
    gadgetron_distributed_gadgets_export.h 
    DistributeGadget.h
    DistributeGadget.cpp
    CollectGadget.h
    CollectGadget.cpp
    IsmrmrdAcquisitionDistributeGadget.h
    IsmrmrdAcquisitionDistributeGadget.cpp
    IsmrmrdImageDistributeGadget.h
    IsmrmrdImageDistributeGadget.cpp
)

set_target_properties(gadgetron_distributed PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})                                                                                                                                                                                                      

target_link_libraries(gadgetron_distributed
    gadgetron_gadgetbase
    gadgetron_toolbox_log
    gadgetron_toolbox_cloudbus
    gadgetron_toolbox_gadgettools
    ${ACE_LIBRARIES}
)

install(FILES 
    gadgetron_distributed_gadgets_export.h
    DistributeGadget.h
    CollectGadget.h
    IsmrmrdAcquisitionDistributeGadget.h
    IsmrmrdImageDistributeGadget.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_distributed DESTINATION lib COMPONENT main)
install(FILES config/distributed_default.xml config/distributed_image_default.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
