IF (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_UTIL_GADGETS__)
ENDIF (WIN32)

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core
)

add_library(gadgetron_util SHARED 
    gadgetron_util_gadgets_export.h 
    ParameterRelayGadget.h
    ParameterRelayGadget.cpp
)

set_target_properties(gadgetron_util PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})                                                                                                                                                                                                      

target_link_libraries(gadgetron_util
    gadgetron_gadgetbase
    gadgetron_toolbox_log
    ${ACE_LIBRARIES}
)

install(FILES 
    gadgetron_util_gadgets_export.h
    ParameterRelayGadget.h
    DESTINATION include COMPONENT main)

install(TARGETS gadgetron_util DESTINATION lib COMPONENT main)
