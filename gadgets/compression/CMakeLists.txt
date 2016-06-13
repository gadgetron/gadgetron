find_package(ZFP REQUIRED)
find_package(ISMRMRD REQUIRED)

IF (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_COMPRESSION_GADGETS__)
ENDIF (WIN32)

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/cloudbus
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
    ${ZFP_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
)

add_library(gadgetron_compression SHARED 
    gadgetron_compression_gadgets_export.h 
    GadgetIsmrmrdCompressedReadWrite.h
    GadgetIsmrmrdCompressedReadWrite.cpp
)

set_target_properties(gadgetron_compression PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})                                                                                                                                                                                                      

target_link_libraries(gadgetron_compression
    gadgetron_gadgetbase
    gadgetron_toolbox_log
    gadgetron_toolbox_gadgettools
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
    ${ZFP_LIBRARIES}
    ${ISMRMRD_LIBRARIES} 
)

install(FILES
    gadgetron_compression_gadgets_export.h  
    GadgetIsmrmrdCompressedReadWrite.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_compression DESTINATION lib COMPONENT main)
install(FILES
        config/default_compression.xml
        config/Generic_Cartesian_Grappa_SNR_Compression.xml
        config/ismrmrd_dump_compressed.xml
        DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
