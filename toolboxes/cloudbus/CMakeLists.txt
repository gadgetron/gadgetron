find_package(ACE REQUIRED)
find_package(Boost REQUIRED)

include_directories(${ACE_INCLUDE_DIR}
                    ${Boost_INCLUDE_DIR} 
                    )

add_library(gadgetron_toolbox_cloudbus SHARED
  CloudBus.cpp
  CloudBus.h
  cloudbus_export.h
)

target_link_libraries(gadgetron_toolbox_cloudbus
		     gadgetron_toolbox_log
                     optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

set_target_properties(gadgetron_toolbox_cloudbus PROPERTIES COMPILE_DEFINITIONS "__BUILD_GADGETRON_CLOUDBUS__")
set_target_properties(gadgetron_toolbox_cloudbus  PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

add_executable(gadgetron_cloudbus cloudbus_main.cpp)
target_link_libraries(gadgetron_cloudbus 
                     gadgetron_toolbox_cloudbus 
		     gadgetron_toolbox_log
                     optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} )

install(TARGETS gadgetron_toolbox_cloudbus DESTINATION lib COMPONENT cloudbus-lib)
install(TARGETS gadgetron_cloudbus DESTINATION bin COMPONENT cloudbus-utils)
install(FILES CloudBus.h cloudbus_export.h DESTINATION include COMPONENT cloudbus-dev)
