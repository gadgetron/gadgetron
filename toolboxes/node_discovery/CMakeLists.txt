if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_NODE_DISCOVERY__)
endif ()

include_directories( ${Boost_INCLUDE_DIR}  )

add_library(gadgetron_toolbox_node_discovery SHARED node_discovery_export.h node_discovery.h node_discovery.cpp)
target_link_libraries(gadgetron_toolbox_node_discovery gadgetron_toolbox_log ${Boost_LIBRARIES})
set_target_properties(gadgetron_toolbox_node_discovery PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS gadgetron_toolbox_node_discovery DESTINATION lib COMPONENT main)
install(FILES node_discovery.h node_discovery_export.h DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

