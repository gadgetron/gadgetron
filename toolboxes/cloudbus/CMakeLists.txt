find_package(ACE REQUIRED)
find_package(Boost REQUIRED)

if (WIN32)
  add_definitions(-DCROW_MSVC_WORKAROUND)
endif()

include_directories(${ACE_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIR}
		    ${CMAKE_SOURCE_DIR}/toolboxes/rest
		    ${CMAKE_SOURCE_DIR}/toolboxes/node_discovery
                    )

add_library(gadgetron_toolbox_cloudbus SHARED
  CloudBus.cpp
  CloudBus.h
  cloudbus_io.h 
  cloudbus_io.cpp 
  cloudbus_export.h
)

target_link_libraries(gadgetron_toolbox_cloudbus
                      gadgetron_toolbox_log
                      gadgetron_toolbox_node_discovery  
                     optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

if(WIN32)
    target_link_libraries(gadgetron_toolbox_cloudbus bcrypt)
endif()

set_target_properties(gadgetron_toolbox_cloudbus PROPERTIES COMPILE_DEFINITIONS "__BUILD_GADGETRON_CLOUDBUS__")
set_target_properties(gadgetron_toolbox_cloudbus  PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

add_executable(gadgetron_cloudbus cloudbus_main.cpp)
target_link_libraries(gadgetron_cloudbus 
                     gadgetron_toolbox_cloudbus 
		     gadgetron_toolbox_log
                     optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} )

add_executable(gadgetron_cloudbus_relay cloudbus_relay.cpp)

target_link_libraries(gadgetron_cloudbus_relay 
                     gadgetron_toolbox_cloudbus 
		     gadgetron_toolbox_log
		     gadgetron_toolbox_rest
                     optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} )

install(TARGETS gadgetron_toolbox_cloudbus DESTINATION lib COMPONENT main)
install(TARGETS gadgetron_cloudbus DESTINATION bin COMPONENT main)
install(TARGETS gadgetron_cloudbus_relay DESTINATION bin COMPONENT main)

install(FILES 
  CloudBus.h
  cloudbus_io.h 
  cloudbus_export.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
