if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_DISTRIBUTED_GADGETS__)
endif ()


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
)

install(FILES 
    gadgetron_distributed_gadgets_export.h
    DistributeGadget.h
    CollectGadget.h
    IsmrmrdAcquisitionDistributeGadget.h
    IsmrmrdImageDistributeGadget.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_distributed DESTINATION lib COMPONENT main)
install(FILES config/distributed_default_old.xml config/distributed_image_default.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
